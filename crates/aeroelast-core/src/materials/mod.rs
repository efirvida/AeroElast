pub mod composite;
pub mod failure;
pub mod isotropic;
pub mod laminate;
pub mod orthotropic;

use nalgebra::{Matrix2, Matrix3};

/// Element family classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementFamily {
    Shell = 2,
    Plane = 3,
}

/// Constitutive matrices for a shell element.
///
/// All matrices are in local element coordinates.
#[derive(Clone)]
pub struct ShellConstitutive {
    /// Membrane stiffness (3×3): force-strain [N·h / (1-ν²)] form
    pub cm: Matrix3<f64>,
    /// Coupling stiffness (3×3): membrane-bending coupling [N] form
    /// N = A·ε + B·κ, M = B·ε + D·κ
    pub cb_coupling: Matrix3<f64>,
    /// Bending stiffness (3×3): moment-curvature [N·h³/12 / (1-ν²)] form
    pub cb: Matrix3<f64>,
    /// Transverse shear stiffness (2×2): force-shear strain [k·G·h] form
    pub cs: Matrix2<f64>,
    /// Raw membrane (stress-strain, no thickness): σ = C_raw · ε
    pub cm_raw: Matrix3<f64>,
}

impl ShellConstitutive {
    /// The **uncorrected** transverse-shear stiffness `G·h`.
    ///
    /// Ko, Bathe & Zhang (2025), *Computers & Structures* 308:107622, p. 410: the
    /// element formulation "does not include any numerical factor", so it consumes
    /// the uncorrected transverse-shear stiffness. `cs` is in `k·G·h` form (see the
    /// field doc), but `ShellConstitutive` does not carry `k`; the caller supplies
    /// the scalar `applied_k` that the material model actually applied when it built
    /// `cs`. `applied_k = 1.0` means "the material applied no scalar factor" — e.g. a
    /// multi-ply laminate's energy-equivalent section stiffness — and `cs` then passes
    /// through verbatim. A non-positive `applied_k` is treated the same way, so the
    /// method never divides by zero.
    ///
    /// Additive: `cs` and every existing consumer are untouched.
    pub fn transverse_shear_uncorrected(&self, applied_k: f64) -> Matrix2<f64> {
        if applied_k > 0.0 { self.cs / applied_k } else { self.cs }
    }
}

/// Material trait for shell elements.
///
/// Implementors provide constitutive matrices and physical properties.
pub trait Material: Send + Sync {
    /// Compute all constitutive matrices for a given thickness and shear correction factor.
    fn constitutive(&self, thickness: f64, shear_correction: f64) -> ShellConstitutive;

    /// Density (kg/m³)
    fn density(&self) -> f64;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::materials::isotropic::IsotropicMaterial;
    use crate::materials::laminate::{Laminate, Ply};
    use crate::materials::orthotropic::OrthotropicMaterial;

    const E: f64 = 70.0e9;
    const NU: f64 = 0.3;
    /// Plane-stress shear modulus `G = E / (2 (1 + ν))`.
    const G: f64 = E / (2.0 * (1.0 + NU));
    const K: f64 = 5.0 / 6.0;

    /// A single-ply laminate whose ply is isotropic-equivalent, so its
    /// transverse-shear block is exactly the scalar-corrected `k·G·h`.
    fn single_isotropic_ply(h: f64, k: f64) -> Laminate {
        let mat = OrthotropicMaterial::new(E, E, E, G, G, G, NU, NU, NU, 2700.0);
        Laminate::new(vec![Ply::new(mat, h, 0.0)], k).unwrap()
    }

    /// A three-ply laminate, i.e. the multi-ply energy-equivalence branch of
    /// `Laminate::compute_shear_stiffness`, which applies no scalar factor.
    fn multi_ply(h: f64, k: f64) -> Laminate {
        let mat = OrthotropicMaterial::new(E, E, E, G, G, G, NU, NU, NU, 2700.0);
        Laminate::new(
            vec![Ply::new(mat, h, 0.0), Ply::new(mat, h, 90.0), Ply::new(mat, h, 0.0)],
            k,
        )
        .unwrap()
    }

    /// A single isotropic ply really is `k·G·h`, so the accessor must hand back
    /// the uncorrected `G·h` once the applied factor is removed.
    #[test]
    fn test_transverse_shear_uncorrected_single_isotropic_ply_equals_g_h() {
        let h = 0.01;
        let lam = single_isotropic_ply(h, K);
        let shell = lam.to_shell_constitutive();

        // The single-ply branch builds `cs = k·G·h`.
        let expected_cs = Matrix2::new(K * G * h, 0.0, 0.0, K * G * h);
        let cs_err = (shell.cs - expected_cs).norm() / expected_cs.norm();
        assert!(cs_err <= 1e-12, "single-ply cs should be k·G·h, rel err {cs_err}");

        // The reported applied factor is the laminate's scalar.
        let applied_k = lam.applied_shear_correction_factor();
        assert!((applied_k - K).abs() <= 1e-15, "applied k = {applied_k}, expected {K}");

        // Removing it recovers the uncorrected `G·h`.
        let uncorrected = shell.transverse_shear_uncorrected(applied_k);
        let expected = Matrix2::new(G * h, 0.0, 0.0, G * h);
        let err = (uncorrected - expected).norm() / expected.norm();
        assert!(err <= 1e-12, "uncorrected should be G·h, rel err {err}");
    }

    /// The isotropic constitutive is also `k·G·h`, so the same removal applies.
    #[test]
    fn test_transverse_shear_uncorrected_isotropic_constitutive_removes_k() {
        let h = 0.01;
        let shell = IsotropicMaterial::new(E, NU, 2700.0).constitutive(h, K);
        let uncorrected = shell.transverse_shear_uncorrected(K);
        let expected = Matrix2::new(G * h, 0.0, 0.0, G * h);
        let err = (uncorrected - expected).norm() / expected.norm();
        assert!(err <= 1e-12, "uncorrected should be G·h, rel err {err}");
    }

    /// A multi-ply laminate's `cs` is the energy-equivalent section stiffness
    /// with no scalar factor, so the accessor must pass it through verbatim.
    #[test]
    fn test_transverse_shear_uncorrected_multi_ply_is_cs_unchanged() {
        let h = 0.01;
        let lam = multi_ply(h, K);
        let shell = lam.to_shell_constitutive();

        assert_eq!(lam.applied_shear_correction_factor(), 1.0);

        let uncorrected = shell.transverse_shear_uncorrected(lam.applied_shear_correction_factor());
        let err = (uncorrected - shell.cs).norm() / shell.cs.norm();
        assert!(err <= 1e-12, "multi-ply uncorrected must equal cs, rel err {err}");
    }

    /// The discriminator: the rejected mechanism (c) — divide the multi-ply
    /// `cs` by the laminate's scalar `shear_correction_factor` — must be
    /// separated from the correct value by more than `1e-3` relative, so the
    /// test cannot pass vacuously.
    #[test]
    fn test_transverse_shear_uncorrected_discriminates_naive_multi_ply_division() {
        let h = 0.01;
        let lam = multi_ply(h, K);
        let shell = lam.to_shell_constitutive();
        let correct = shell.transverse_shear_uncorrected(lam.applied_shear_correction_factor());

        let wrong = shell.cs / lam.shear_correction_factor;
        let rel = (wrong - correct).norm() / correct.norm();
        assert!(rel > 1e-3, "naive division must be separated, rel = {rel}");
    }

    /// A non-positive factor means "the material applied none": `cs` passes
    /// through unchanged, with no division by zero.
    #[test]
    fn test_transverse_shear_uncorrected_nonpositive_factor_returns_cs() {
        let shell = IsotropicMaterial::new(E, NU, 2700.0).constitutive(0.01, K);
        assert_eq!(shell.transverse_shear_uncorrected(0.0), shell.cs);
        assert_eq!(shell.transverse_shear_uncorrected(-1.0), shell.cs);
    }
}
