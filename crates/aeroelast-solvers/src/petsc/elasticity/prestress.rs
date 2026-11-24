//! Centrifugal static-solve prestress and its geometric stiffness `K_G`.
//!
//! One implementation of the physics that the rotating-blade solvers share:
//!
//! 1. nodal centrifugal load `f_cf` (lumped per element, mass per area),
//! 2. static solve `K u = f_cf` on the free DOFs (linear `K`, clamped-root
//!    elimination), reduced `K` factorized once and cached,
//! 3. membrane resultants `(A eps + B kappa)/h` recovered at the element
//!    centroids from `u` (`compute_prestress_field`),
//! 4. `K_G = assemble_geometric_k(sigma)`.
//!
//! The prestress is always solved at the reference coordinates the assembler
//! holds when the call starts (`X0`) with the linear `K` ("Option A
//! discipline").  With `deformed_coords = Some(X0 + u)` only the final
//! `assemble_geometric_k` runs on those coordinates; `X0` is restored before
//! returning.
//!
//! The old local estimate `sigma = rho w^2 r l_char` (`assemble_centrifugal_k`)
//! was two to three orders of magnitude too small and is gone, with no
//! fallback.  Python `MeshAssembler.assemble_geometric_stiffness` reaches this
//! code through the PyO3 binding, so there is a single K_G for centrifugal
//! prestress (L33, T17).

use aeroelast_core::assembly::assembler::MeshAssembler;

use super::super::assembler::{assemble_seq_aij, create_vec, ensure_initialized};
use super::super::infra::ffi::{self, INSERT_VALUES};
use super::super::infra::handles::PetscKsp;
use super::super::infra::mat::{check, PetscError, PetscMat};
use super::super::infra::vec::PetscVec;

const KSPPREONLY: &std::ffi::CStr =
    unsafe { std::ffi::CStr::from_bytes_with_nul_unchecked(b"preonly\0") };
const PCLU: &std::ffi::CStr =
    unsafe { std::ffi::CStr::from_bytes_with_nul_unchecked(b"lu\0") };

/// Cached factorization of the reduced linear `K` plus the K_G recipe.
pub struct CentrifugalPrestress {
    // Field order matters only for readability: PETSc reference-counts the Mat.
    ksp: PetscKsp,
    _k_red: PetscMat,
    /// Free global DOFs; the reduced index of `free_dofs[i]` is `i`.
    free_dofs: Vec<usize>,
    n_full: usize,
}

fn reduced_k(asm: &MeshAssembler, free_dofs: &[usize]) -> Result<PetscMat, PetscError> {
    let n_full = asm.dofs_count;
    let mut full_to_red = vec![-1i32; n_full];
    for (i, &d) in free_dofs.iter().enumerate() {
        assert!(d < n_full, "free DOF {d} outside 0..{n_full}");
        full_to_red[d] = i as i32;
    }
    let (rows, cols, vals) = asm.assemble_k();
    let mut r_red = Vec::with_capacity(vals.len());
    let mut c_red = Vec::with_capacity(vals.len());
    let mut v_red = Vec::with_capacity(vals.len());
    for ((&r, &c), &v) in rows.iter().zip(cols.iter()).zip(vals.iter()) {
        let (ri, ci) = (full_to_red[r as usize], full_to_red[c as usize]);
        if ri >= 0 && ci >= 0 {
            r_red.push(ri);
            c_red.push(ci);
            v_red.push(v);
        }
    }
    assemble_seq_aij(&r_red, &c_red, &v_red, free_dofs.len())
}

fn build_vec(values: &[f64]) -> Result<PetscVec, PetscError> {
    let n = values.len();
    let v = create_vec(n)?;
    let indices: Vec<i32> = (0..n as i32).collect();
    unsafe {
        check(
            ffi::VecSetValues(v.as_raw(), n as i32, indices.as_ptr(), values.as_ptr(), INSERT_VALUES),
            "VecSetValues",
        )?;
        check(ffi::VecAssemblyBegin(v.as_raw()), "VecAssemblyBegin")?;
        check(ffi::VecAssemblyEnd(v.as_raw()), "VecAssemblyEnd")?;
    }
    Ok(v)
}

fn create_ksp(k_red: &PetscMat) -> Result<PetscKsp, PetscError> {
    ensure_initialized()?;
    unsafe {
        let mut raw: ffi::KSP = std::ptr::null_mut();
        check(ffi::KSPCreate(ffi::petsc_comm_self(), &mut raw), "KSPCreate")?;
        let ksp = PetscKsp::from_raw(raw);
        check(ffi::KSPSetType(ksp.as_raw(), KSPPREONLY.as_ptr()), "KSPSetType(preonly)")?;
        let mut pc: ffi::PC = std::ptr::null_mut();
        check(ffi::KSPGetPC(ksp.as_raw(), &mut pc), "KSPGetPC")?;
        check(ffi::PCSetType(pc, PCLU.as_ptr()), "PCSetType(lu)")?;
        check(
            ffi::KSPSetOperators(ksp.as_raw(), k_red.as_raw(), k_red.as_raw()),
            "KSPSetOperators",
        )?;
        check(ffi::KSPSetFromOptions(ksp.as_raw()), "KSPSetFromOptions")?;
        Ok(ksp)
    }
}

impl CentrifugalPrestress {
    /// Factorization of the linear `K` of `asm` (at its current coordinates)
    /// on `free_dofs`.  PETSc factorizes lazily on the first solve.
    pub fn new(asm: &MeshAssembler, free_dofs: &[usize]) -> Result<Self, PetscError> {
        let k_red = reduced_k(asm, free_dofs)?;
        let ksp = create_ksp(&k_red)?;
        Ok(Self {
            ksp,
            _k_red: k_red,
            free_dofs: free_dofs.to_vec(),
            n_full: asm.dofs_count,
        })
    }

    /// Re-read the linear `K` from `asm` (its coordinates changed, e.g. the
    /// inertial solver's `X_rotated(theta)`); refactorizes on the next solve.
    pub fn refresh_operator(&mut self, asm: &MeshAssembler) -> Result<(), PetscError> {
        assert_eq!(asm.dofs_count, self.n_full, "assembler DOF count changed");
        let k_red = reduced_k(asm, &self.free_dofs)?;
        let ksp = create_ksp(&k_red)?;
        self.ksp = ksp;
        self._k_red = k_red;
        Ok(())
    }

    /// Static prestress displacement `K u = f_cf` (full DOF vector, zero on
    /// the constrained DOFs).  `mass_per_area` is `asm.mass_per_area_per_elem()`.
    pub fn displacement(
        &self,
        asm: &MeshAssembler,
        omega: f64,
        axis: [f64; 3],
        center: [f64; 3],
    ) -> Result<Vec<f64>, PetscError> {
        let mpa = asm.mass_per_area_per_elem();
        let f_full = asm.centrifugal_load(omega, axis, center, &mpa);
        let f_red: Vec<f64> = self.free_dofs.iter().map(|&d| f_full[d]).collect();
        let rhs = build_vec(&f_red)?;
        let sol = create_vec(self.free_dofs.len())?;
        unsafe {
            check(ffi::KSPSolve(self.ksp.as_raw(), rhs.as_raw(), sol.as_raw()), "KSPSolve")?;
            let mut reason: i32 = 0;
            check(
                ffi::KSPGetConvergedReason(self.ksp.as_raw(), &mut reason),
                "KSPGetConvergedReason",
            )?;
            if reason <= 0 {
                return Err(PetscError {
                    code: reason,
                    context: "centrifugal prestress solve K u = f_cf did not converge",
                });
            }
        }
        let u_red = sol.to_vec()?;
        let mut u_full = vec![0.0; self.n_full];
        for (i, &d) in self.free_dofs.iter().enumerate() {
            u_full[d] = u_red[i];
        }
        Ok(u_full)
    }

    /// Centrifugal `K_G` COO triplets (full DOF numbering).
    ///
    /// The prestress is solved and recovered on the coordinates `asm` holds
    /// now (`X0`).  `deformed_coords = Some(X0 + u)` makes only the final
    /// `assemble_geometric_k` use those coordinates; `X0` is restored.
    pub fn geometric_k(
        &self,
        asm: &mut MeshAssembler,
        omega: f64,
        axis: [f64; 3],
        center: [f64; 3],
        deformed_coords: Option<&[f64]>,
    ) -> Result<(Vec<i64>, Vec<i64>, Vec<f64>), PetscError> {
        let u = self.displacement(asm, omega, axis, center)?;
        let sigma = asm.compute_prestress_field(&u);
        match deformed_coords {
            None => Ok(asm.assemble_geometric_k(&sigma)),
            Some(coords) => {
                let x0 = asm.topology.node_coords.clone();
                asm.update_node_coordinates(coords);
                let out = asm.assemble_geometric_k(&sigma);
                asm.update_node_coordinates(&x0);
                Ok(out)
            }
        }
    }
}

/// Shared test fixture: clamped isotropic plate in the x-y plane, hub offset.
#[cfg(test)]
pub(crate) mod test_support {
    use aeroelast_core::assembly::assembler::{MaterialSpec, MeshAssembler};
    use aeroelast_core::assembly::topology::{ElemType, MeshTopology};

    /// `nx` x `ny` MITC4 quads, nodes at `x = hub + i`, `y = j`; returns the
    /// assembler and the free DOFs (clamp at `x == hub`, all 6 DOFs).
    pub(crate) fn hub_plate(nx: usize, ny: usize, hub: f64) -> (MeshAssembler, Vec<usize>) {
        let nid = |i: usize, j: usize| j * (nx + 1) + i;
        let mut coords = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                coords.extend_from_slice(&[hub + i as f64, j as f64, 0.0]);
            }
        }
        let mut conn = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                conn.push(vec![nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]);
            }
        }
        let n_el = conn.len();
        let topo = MeshTopology::new(coords, conn, vec![ElemType::Mitc4; n_el]);
        let mat = MaterialSpec::Isotropic {
            e: 2.0e11,
            nu: 0.3,
            rho: 7800.0,
            thickness: 0.02,
            shear_correction: 5.0 / 6.0,
            drilling_scale: 1.0,
        };
        let asm = MeshAssembler::new(topo, vec![mat; n_el]);
        let mut free = Vec::new();
        for j in 0..=ny {
            for i in 1..=nx {
                for d in 0..6 {
                    free.push(6 * nid(i, j) + d);
                }
            }
        }
        free.sort_unstable();
        (asm, free)
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::hub_plate;
    use super::*;

    const AXIS: [f64; 3] = [0.0, 0.0, 1.0];
    const CENTER: [f64; 3] = [0.0; 3];

    fn dense(rows: &[i64], cols: &[i64], vals: &[f64], n: usize) -> Vec<f64> {
        let mut m = vec![0.0; n * n];
        for ((&r, &c), &v) in rows.iter().zip(cols).zip(vals) {
            m[r as usize * n + c as usize] += v;
        }
        m
    }

    #[test]
    fn mass_per_area_is_rho_times_thickness_for_isotropic() {
        let (asm, _) = hub_plate(2, 1, 1.0);
        for m in asm.mass_per_area_per_elem() {
            assert!((m - 7800.0 * 0.02).abs() < 1e-12);
        }
    }

    #[test]
    fn displacement_is_radially_tensile_and_equilibrated() {
        let (asm, free) = hub_plate(4, 2, 1.0);
        let cp = CentrifugalPrestress::new(&asm, &free).unwrap();
        let u = cp.displacement(&asm, 12.0, AXIS, CENTER).unwrap();
        // K u = f on the free DOFs (equilibrium of the static solve).
        let mpa = asm.mass_per_area_per_elem();
        let f = asm.centrifugal_load(12.0, AXIS, CENTER, &mpa);
        let (r, c, v) = asm.assemble_k();
        let mut ku = vec![0.0; asm.dofs_count];
        for ((&r, &c), &v) in r.iter().zip(&c).zip(&v) {
            ku[r as usize] += v * u[c as usize];
        }
        let fmax = f.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        for &d in &free {
            assert!((ku[d] - f[d]).abs() < 1e-8 * fmax, "dof {d}");
        }
        // Clamped DOFs untouched; the y = 0 edge tip (node nx = 4) moves outward (+x).
        assert_eq!(u[0], 0.0);
        assert!(u[6 * 4] > 0.0, "tip u_x = {:e}", u[6 * 4]);
    }

    #[test]
    fn geometric_k_nonzero_symmetric_and_restores_coords() {
        let (mut asm, free) = hub_plate(4, 2, 1.0);
        let cp = CentrifugalPrestress::new(&asm, &free).unwrap();
        let x0 = asm.topology.node_coords.clone();
        let (r, c, v) = cp.geometric_k(&mut asm, 12.0, AXIS, CENTER, None).unwrap();
        let n = asm.dofs_count;
        let kg = dense(&r, &c, &v, n);
        assert!(kg.iter().any(|x| x.abs() > 1.0), "K_G must be non-trivial");
        for i in 0..n {
            for j in 0..n {
                let (a, b) = (kg[i * n + j], kg[j * n + i]);
                assert!((a - b).abs() <= 1e-9 * (a.abs() + b.abs() + 1.0));
            }
        }
        // Deformed-coords variant: prestress still at X0, final assembly on X0+u.
        let u: Vec<f64> = (0..n).map(|i| if i % 6 == 0 { 1e-3 } else { 0.0 }).collect();
        let xd = aeroelast_core_compose(&x0, &u);
        let (r2, c2, v2) = cp.geometric_k(&mut asm, 12.0, AXIS, CENTER, Some(&xd)).unwrap();
        assert_eq!(asm.topology.node_coords, x0, "X0 must be restored");
        assert_eq!(r2.len(), r.len());
        assert_eq!(c2.len(), c.len());
        assert!(v2.iter().zip(&v).any(|(a, b)| a != b), "X0+u must change K_G");
    }

    fn aeroelast_core_compose(x0: &[f64], u: &[f64]) -> Vec<f64> {
        x0.chunks(3)
            .enumerate()
            .flat_map(|(i, p)| (0..3).map(move |d| p[d] + u[6 * i + d]))
            .collect()
    }

    #[test]
    fn refresh_operator_tracks_new_coordinates() {
        let (mut asm, free) = hub_plate(2, 1, 1.0);
        let mut cp = CentrifugalPrestress::new(&asm, &free).unwrap();
        let u_a = cp.displacement(&asm, 5.0, AXIS, CENTER).unwrap();
        // Stretch the mesh 2x in x; the operator must be refreshed to stay consistent.
        let stretched: Vec<f64> = asm
            .topology
            .node_coords
            .chunks(3)
            .flat_map(|p| [2.0 * p[0], p[1], p[2]])
            .collect();
        asm.update_node_coordinates(&stretched);
        cp.refresh_operator(&asm).unwrap();
        let fresh = CentrifugalPrestress::new(&asm, &free).unwrap();
        let u_b = cp.displacement(&asm, 5.0, AXIS, CENTER).unwrap();
        let u_c = fresh.displacement(&asm, 5.0, AXIS, CENTER).unwrap();
        assert_eq!(u_b, u_c);
        assert_ne!(u_a, u_b);
    }
}
