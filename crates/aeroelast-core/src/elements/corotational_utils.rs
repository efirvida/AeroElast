use nalgebra::Matrix3;

/// Element-level orthonormal frame used by corotational shell primitives.
#[derive(Clone, Copy, Debug)]
pub struct LocalFrame3 {
    pub e1: nalgebra::Vector3<f64>,
    pub e2: nalgebra::Vector3<f64>,
    pub e3: nalgebra::Vector3<f64>,
}

/// Polar decomposition: F = R · U via SVD, with det(R)=+1 correction.
///
/// Input is the displacement gradient `H`, with `F = I + H`.
pub fn polar_decomposition(h: &Matrix3<f64>) -> (Matrix3<f64>, Matrix3<f64>) {
    let f = Matrix3::identity() + h;

    let svd = f.svd(true, true);
    let u_mat = svd.u.expect("SVD U failed");
    let vt = svd.v_t.expect("SVD Vt failed");

    let mut r = u_mat * vt;

    if r.determinant() < 0.0 {
        let mut min_idx = 0usize;
        let mut min_sigma = svd.singular_values[0];
        for i in 1..3 {
            if svd.singular_values[i] < min_sigma {
                min_sigma = svd.singular_values[i];
                min_idx = i;
            }
        }

        let mut u_corr = u_mat;
        u_corr.set_column(min_idx, &(-u_mat.column(min_idx)));
        r = u_corr * vt;
    }

    let u_stretch = r.transpose() * f;
    (r, u_stretch)
}
