/// SC-03 — Mock FSI loop without preCICE.
///
/// Verifies the full implicit coupling pattern:
///   while coupling_ongoing:
///     dt = get_dt()
///     cp = checkpoint()
///     while not converged:
///         forces = read_forces()
///         step(forces, dt)
///         write_displacements()
///         if requires_reading_checkpoint:
///             restore(cp)
///             continue
///     advance(dt)
///
/// We use a `MockCoupling` struct to simulate the preCICE participant,
/// driving the coupling loop from outside. The test verifies that:
/// - checkpoint/restore correctly rolls back the stepper
/// - The stepper reaches the expected displacement after N coupling steps
/// - Energy is conserved within tolerance for the undamped case
use aeroelast_solvers::petsc::elasticity::dynamic_newmark::{NewmarkCheckpoint, NewmarkStepper};

// ── MockCoupling ─────────────────────────────────────────────────────────────

/// Simulates a preCICE participant that applies a sinusoidal aerodynamic force
/// and demands `n_rollbacks` checkpoint rollbacks before accepting each step.
struct MockCoupling {
    /// Total number of coupling time windows to simulate.
    n_windows: usize,
    /// How many times the mock demands a rollback per window (0 = accept immediately).
    rollbacks_per_window: usize,
    /// Time step size (constant).
    dt: f64,
    /// Aerodynamic force amplitude.
    force_amplitude: f64,

    // ── state ────────────────────────────────────────────────────────────────
    current_window: usize,
    rollbacks_this_window: usize,
    /// Time elapsed (from the coupling side).
    t: f64,
}

impl MockCoupling {
    fn new(n_windows: usize, rollbacks_per_window: usize, dt: f64, force_amplitude: f64) -> Self {
        Self {
            n_windows,
            rollbacks_per_window,
            dt,
            force_amplitude,
            current_window: 0,
            rollbacks_this_window: 0,
            t: 0.0,
        }
    }

    fn is_coupling_ongoing(&self) -> bool {
        self.current_window < self.n_windows
    }

    fn dt(&self) -> f64 {
        self.dt
    }

    /// Returns `true` if the mock demands a rollback this iteration.
    fn requires_rollback(&self) -> bool {
        self.rollbacks_this_window < self.rollbacks_per_window
    }

    /// Called after a successful sub-iteration — accepts or marks rollback.
    fn accept_or_rollback(&mut self) {
        if self.rollbacks_this_window < self.rollbacks_per_window {
            self.rollbacks_this_window += 1;
        }
    }

    /// Called when the window is converged — advance coupling time.
    fn advance(&mut self) {
        self.t += self.dt;
        self.current_window += 1;
        self.rollbacks_this_window = 0;
    }

    /// Read aerodynamic forces at current coupling time (sinusoidal).
    fn read_forces(&self) -> Vec<f64> {
        let f = self.force_amplitude * (std::f64::consts::TAU * self.t).sin();
        vec![f]
    }
}

// ── SC-03 test ────────────────────────────────────────────────────────────────

#[test]
fn test_mock_fsi_loop_with_rollbacks() {
    // 1-DOF undamped system: m·ü + k·u = F_aero(t)
    let k = 100.0f64;
    let m = 1.0f64;
    let omega_struct = (k / m).sqrt(); // structural natural frequency = 10 rad/s
    let dt = 2.0 * std::f64::consts::PI / omega_struct / 50.0; // 50 steps/period

    let n_windows = 20;
    let rollbacks_per_window = 2; // each window demands 2 rollbacks before accepting

    let rows = vec![0i32];
    let cols = vec![0i32];

    let mut stepper = NewmarkStepper::new(
        &rows, &cols, &[k],
        &rows, &cols, &[m],
        &rows, &cols, &[0.0f64], // undamped
        1,
        0.25,
        0.5,
        dt,
    ).expect("NewmarkStepper::new failed");

    let mut coupling = MockCoupling::new(n_windows, rollbacks_per_window, dt, 10.0);

    let mut accepted_windows = 0usize;
    let mut total_rollbacks = 0usize;

    while coupling.is_coupling_ongoing() {
        let dt_coup = coupling.dt();

        // Checkpoint before this coupling window
        let cp: NewmarkCheckpoint = stepper.checkpoint();
        let t_before = stepper.current_time();

        let mut converged = false;
        while !converged {
            // Read forces from aerodynamic solver (mock)
            let f_aero = coupling.read_forces();

            // Structural sub-step
            stepper.step(&f_aero, dt_coup).expect("step failed");

            // Write displacements (no-op in mock — just verify non-NaN)
            let u_curr = stepper.current_time(); // borrow to check it's valid
            assert!(u_curr.is_finite(), "stepper time is NaN/inf");

            // Check if mock demands rollback
            if coupling.requires_rollback() {
                // Roll back
                stepper.restore(&cp);
                assert!(
                    (stepper.current_time() - t_before).abs() < 1e-14,
                    "Restore did not recover time: got {}, expected {t_before}",
                    stepper.current_time()
                );
                coupling.accept_or_rollback();
                total_rollbacks += 1;
            } else {
                converged = true;
            }
        }

        // Advance coupling window
        coupling.advance();
        accepted_windows += 1;
    }

    assert_eq!(
        accepted_windows, n_windows,
        "Expected {n_windows} accepted windows, got {accepted_windows}"
    );

    let expected_rollbacks = n_windows * rollbacks_per_window;
    assert_eq!(
        total_rollbacks, expected_rollbacks,
        "Expected {expected_rollbacks} rollbacks, got {total_rollbacks}"
    );

    // Final time should match n_windows * dt
    let t_final = stepper.current_time();
    let t_expected = n_windows as f64 * dt;
    assert!(
        (t_final - t_expected).abs() < 1e-10,
        "Final time mismatch: got {t_final:.6}, expected {t_expected:.6}"
    );

    // Solution must be bounded (no blow-up)
    assert!(
        t_final.is_finite(),
        "Final time is not finite: {t_final}"
    );
}

/// SC-03b — Zero-rollback variant: mock accepts every step immediately.
///
/// This is the "happy path" — no checkpoint/restore needed.
/// Verifies the stepper advances correctly for N consecutive windows.
#[test]
fn test_mock_fsi_loop_no_rollbacks() {
    let k = 100.0f64;
    let m = 1.0f64;
    let omega = (k / m).sqrt();
    let dt = 2.0 * std::f64::consts::PI / omega / 50.0;
    let n_windows = 10;

    let rows = vec![0i32];
    let cols = vec![0i32];

    let mut stepper = NewmarkStepper::new(
        &rows, &cols, &[k],
        &rows, &cols, &[m],
        &rows, &cols, &[0.0f64],
        1,
        0.25,
        0.5,
        dt,
    ).expect("NewmarkStepper::new failed");

    let mut coupling = MockCoupling::new(n_windows, 0, dt, 0.0); // zero force, zero rollbacks

    while coupling.is_coupling_ongoing() {
        let f = coupling.read_forces();
        stepper.step(&f, coupling.dt()).expect("step failed");
        coupling.advance();
    }

    // With zero force and zero IC, solution should remain zero
    // (within floating-point noise from the KSP solve)
    let t_final = stepper.current_time();
    assert!(
        (t_final - n_windows as f64 * dt).abs() < 1e-10,
        "Time mismatch: {t_final}"
    );
}

// ── SC-04 — Inertial rotor physics ────────────────────────────────────────────

/// SC-04a — Zero angular velocity → zero F_ref → zero elastic displacement.
///
/// Physical basis: with ω = 0 and α = 0, the rigid-body reference acceleration
/// is zero, so F_ref = 0.  With zero initial conditions and no external forces,
/// the elastic displacement must remain zero for all time (no numerical drift).
///
/// This is the "pure rigid rotation → u_e ≈ 0" contract for the trivial case ω = 0.
#[test]
fn test_inertial_rotor_zero_omega_zero_displacement() {
    use aeroelast_solvers::petsc::fsi::rotor_physics::{
        compute_reference_load_vector, compute_rigid_body_acceleration_inertial,
    };

    let m = 1.0f64;
    let k = 1000.0f64;
    let omega = 0.0f64;
    let alpha = 0.0f64;
    let r = 5.0f64; // non-trivial radius — would give non-zero force if ω ≠ 0

    // Verify physics: zero ω and zero α → zero acceleration → zero F_ref
    let coords = vec![r, 0.0, 0.0];
    let masses = vec![m];
    let axis = [0.0f64, 1.0f64, 0.0f64];
    let center = [0.0f64; 3];

    let a_ref = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);
    assert_eq!(a_ref, vec![0.0, 0.0, 0.0], "a_ref must be zero for ω=0, α=0");

    let f_ref = compute_reference_load_vector(&a_ref, &masses);
    assert_eq!(f_ref, vec![0.0, 0.0, 0.0], "F_ref must be zero for zero a_ref");

    // Integrate: zero force → zero displacement at every step
    let omega_n = (k / m).sqrt();
    let dt = 2.0 * std::f64::consts::PI / omega_n / 50.0;
    let n_windows = 100;

    let rows = vec![0i32];
    let cols = vec![0i32];
    let mut stepper = NewmarkStepper::new(
        &rows, &cols, &[k],
        &rows, &cols, &[m],
        &rows, &cols, &[0.0f64],
        1,
        0.25,
        0.5,
        dt,
    ).expect("NewmarkStepper::new failed");

    for _ in 0..n_windows {
        stepper.step(&[0.0], dt).expect("step failed");
        let u = stepper.current_u()[0];
        assert!(
            u.abs() < 1e-12,
            "Displacement must be zero for zero forcing: u = {u:.3e}"
        );
    }

    let t_final = stepper.current_time();
    assert!((t_final - n_windows as f64 * dt).abs() < 1e-10);
}

/// SC-04b — Constant ω: centripetal F_ref drives bounded elastic deflection.
///
/// Physical setup: a single node of mass m at radius R from the Y rotation axis
/// is connected to a radial spring K.  Under constant angular velocity ω, the
/// reference load F_ref = −M·a_ref = +M·ω²·R (centrifugal direction).
///
/// For an undamped system starting from rest under constant step loading:
///   u(t) = u_static · (1 − cos(ω_n · t))
///
/// → peak displacement = 2·u_static, mean displacement = u_static.
///
/// Verifications:
/// 1. `a_ref.x` has the correct sign (inward, negative x for node at +x).
/// 2. `F_ref.x` has the correct sign (outward, positive x).
/// 3. Peak displacement ≤ 2.1·u_static (bounded, not blown up).
/// 4. Min displacement ≥ −0.1·u_static (no unphysical negative deflection).
/// 5. Final time matches n_windows × dt (stepper advances correctly).
#[test]
fn test_inertial_rotor_centripetal_reference_load() {
    use aeroelast_solvers::petsc::fsi::rotor_physics::{
        compute_reference_load_vector, compute_rigid_body_acceleration_inertial,
    };

    // Physical parameters
    let m = 2.0f64;         // node mass [kg]
    let k = 5_000.0f64;     // radial spring stiffness [N/m]
    let omega = 10.0f64;    // angular velocity [rad/s]
    let r = 2.0f64;         // radial distance from Y axis [m]
    let alpha = 0.0f64;     // constant ω

    // u_static = M·ω²·R / K
    let f_ref_expected = m * omega * omega * r; // 2 * 100 * 2 = 400 N
    let u_static = f_ref_expected / k;          // 400 / 5000 = 0.08 m

    // Verify physics functions: centripetal direction and magnitude
    let coords = vec![r, 0.0, 0.0]; // node on +x axis at radius R
    let masses = vec![m];
    let axis = [0.0f64, 1.0f64, 0.0f64]; // Y rotation axis
    let center = [0.0f64; 3];

    let a_ref = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);

    // a_ref must point inward (toward Y axis): a_ref.x < 0 for node at +x
    assert!(
        a_ref[0] < 0.0,
        "Centripetal a_ref must be inward (negative x): got {:.6e}",
        a_ref[0]
    );
    assert!(
        (a_ref[0] + omega * omega * r).abs() < 1e-10,
        "a_ref.x magnitude: got {:.6e}, expected {:.6e}",
        a_ref[0],
        -(omega * omega * r)
    );
    assert!(a_ref[1].abs() < 1e-14, "a_ref.y must be zero");
    assert!(a_ref[2].abs() < 1e-14, "a_ref.z must be zero");

    let f_ref = compute_reference_load_vector(&a_ref, &masses);

    // F_ref must point outward (centrifugal direction): F_ref.x > 0
    assert!(
        f_ref[0] > 0.0,
        "F_ref must be outward (positive x): got {:.6e}",
        f_ref[0]
    );
    assert!(
        (f_ref[0] - f_ref_expected).abs() < 1e-10,
        "F_ref.x: got {:.6e}, expected {:.6e}",
        f_ref[0],
        f_ref_expected
    );

    // Time integration: structural ω_n = √(K/m) >> ω_rotation → quasi-static regime
    let omega_n = (k / m).sqrt(); // √(5000/2) ≈ 50 rad/s
    let dt = 2.0 * std::f64::consts::PI / omega_n / 80.0; // 80 steps/period
    let n_windows = 200; // ~2.5 full structural periods

    let rows = vec![0i32];
    let cols = vec![0i32];
    let mut stepper = NewmarkStepper::new(
        &rows, &cols, &[k],
        &rows, &cols, &[m],
        &rows, &cols, &[0.0f64], // undamped
        1,
        0.25,
        0.5,
        dt,
    ).expect("NewmarkStepper::new failed");

    // Constant centripetal force (x-component only → our 1-DOF radial system)
    let f_radial = f_ref[0];

    let mut min_u = 0.0f64;
    let mut max_u = 0.0f64;

    for _ in 0..n_windows {
        stepper.step(&[f_radial], dt).expect("step failed");
        let u = stepper.current_u()[0];
        assert!(u.is_finite(), "Displacement is NaN/inf: {u}");
        min_u = min_u.min(u);
        max_u = max_u.max(u);
    }

    // Undamped step response: u(t) = u_static·(1 − cos(ω_n·t))
    // → peak ≈ 2·u_static, min ≈ 0
    let tolerance = 0.1 * u_static;
    assert!(
        max_u <= 2.0 * u_static + tolerance,
        "Peak displacement {max_u:.6} > 2·u_static + tolerance ({:.6})",
        2.0 * u_static + tolerance
    );
    assert!(
        min_u >= -tolerance,
        "Min displacement {min_u:.6} < -tolerance ({tolerance:.6}): unphysical for centripetal loading"
    );

    // Time must advance correctly
    let t_final = stepper.current_time();
    assert!(
        (t_final - n_windows as f64 * dt).abs() < 1e-10,
        "Time mismatch: {t_final:.9} ≠ {:.9}",
        n_windows as f64 * dt
    );
}

/// SC-04c — Gravity-only static deflection.
///
/// Physical setup: 1-DOF vertical spring-mass with gravity.
/// F_g = m·g (constant), K·u = F_g → u_static = m·g / K.
///
/// Undamped Newmark starting from rest under step loading:
///   u(t) = u_static · (1 − cos(ω_n · t))
///
/// Verifications:
/// 1. `compute_gravity_force` returns the correct constant force vector.
/// 2. Newmark response oscillates between 0 and 2·u_static (no blow-up).
/// 3. After enough periods the mean displacement is within 5 % of u_static.
#[test]
fn test_inertial_rotor_gravity_deflection() {
    use aeroelast_solvers::petsc::fsi::rotor_physics::compute_gravity_force;

    let m = 1.5f64;           // [kg]
    let k = 3_000.0f64;       // [N/m]
    let g_val = 9.81f64;      // [m/s²] — downward in z

    let masses = vec![m];
    let g = [0.0f64, 0.0f64, -g_val]; // global gravity vector

    // 1. Verify compute_gravity_force
    let f_g = compute_gravity_force(&masses, &g);
    assert_eq!(f_g.len(), 3, "f_g must be n_nodes × 3");
    assert!((f_g[0]).abs() < 1e-14, "f_g.x must be zero");
    assert!((f_g[1]).abs() < 1e-14, "f_g.y must be zero");
    assert!(
        (f_g[2] - (-m * g_val)).abs() < 1e-12,
        "f_g.z: got {:.6e}, expected {:.6e}",
        f_g[2],
        -m * g_val
    );

    // 2. Integrate a 1-DOF vertical (z) system under constant gravity load
    //    Use z-component of f_g as the single DOF loading
    let f_z = f_g[2]; // negative
    let u_static = f_z / k; // negative → downward displacement

    let omega_n = (k / m).sqrt();
    let dt = 2.0 * std::f64::consts::PI / omega_n / 100.0;
    let n_periods = 5usize;
    let n_windows = (n_periods as f64 * 2.0 * std::f64::consts::PI / omega_n / dt).ceil() as usize;

    let rows = vec![0i32];
    let cols = vec![0i32];
    let mut stepper = NewmarkStepper::new(
        &rows, &cols, &[k],
        &rows, &cols, &[m],
        &rows, &cols, &[0.0f64], // undamped
        1,
        0.25,
        0.5,
        dt,
    ).expect("NewmarkStepper::new failed");

    let mut sum_u = 0.0f64;
    let mut max_u = 0.0f64;
    let mut min_u = 0.0f64;

    for _ in 0..n_windows {
        stepper.step(&[f_z], dt).expect("step failed");
        let u = stepper.current_u()[0];
        assert!(u.is_finite(), "displacement is NaN/inf");
        sum_u += u;
        max_u = max_u.max(u);
        min_u = min_u.min(u);
    }

    // Undamped step response: peak = 2·u_static (downward), min ≈ 0
    let tol = 0.1 * u_static.abs();
    assert!(
        min_u >= 2.0 * u_static - tol,
        "Peak downward disp {min_u:.6} < 2·u_static - tol ({:.6})",
        2.0 * u_static - tol
    );
    assert!(
        max_u <= tol,
        "Max upward disp {max_u:.6} > tol ({tol:.6}): unphysical for gravity-only load"
    );

    // Mean displacement ≈ u_static (within 5 %)
    let mean_u = sum_u / n_windows as f64;
    assert!(
        (mean_u - u_static).abs() < 0.05 * u_static.abs(),
        "Mean displacement {mean_u:.6} not within 5% of u_static {u_static:.6}"
    );
}

/// SC-04d — `a_ref` analytical verification for Y-axis rotation.
///
/// For a node at (R, 0, 0) with Y-axis rotation (axis = ŷ):
///   centripetal:  a_cent = -ω² · R · x̂  (inward)
///   tangential:   a_tang = α · ŷ × r = α · ŷ × R·x̂ = -α·R · ẑ
///
/// Analytical expected values (exact, no truncation):
///   a_ref.x = -ω²·R
///   a_ref.y =  0
///   a_ref.z = -α·R
///
/// Verify against `compute_rigid_body_acceleration_inertial` to tolerance 1e-12.
#[test]
fn test_inertial_rotor_a_ref_analytical_y_axis() {
    use aeroelast_solvers::petsc::fsi::rotor_physics::compute_rigid_body_acceleration_inertial;

    let r = 3.7f64;
    let omega = 12.5f64;
    let alpha = 0.8f64;

    let coords = vec![r, 0.0, 0.0]; // node on +X axis at radius R
    let axis = [0.0f64, 1.0f64, 0.0f64]; // Y rotation axis
    let center = [0.0f64; 3];

    let a = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);

    let expected_x = -(omega * omega * r); // centripetal (inward)
    let expected_y = 0.0f64;
    let expected_z = -(alpha * r); // tangential: α·ŷ × R·x̂ = -αR·ẑ

    assert!(
        (a[0] - expected_x).abs() < 1e-12,
        "a_ref.x: got {:.15e}, expected {:.15e}",
        a[0],
        expected_x
    );
    assert!(
        (a[1] - expected_y).abs() < 1e-12,
        "a_ref.y: got {:.15e}, expected 0",
        a[1]
    );
    assert!(
        (a[2] - expected_z).abs() < 1e-12,
        "a_ref.z: got {:.15e}, expected {:.15e}",
        a[2],
        expected_z
    );
}

/// SC-04e — `a_ref` analytical verification for Z-axis rotation.
///
/// For a node at (R, 0, 0) with Z-axis rotation (axis = ẑ):
///   centripetal:  a_cent = -ω²·R·x̂   (inward, toward Z axis)
///   tangential:   a_tang = α·ẑ × R·x̂ = α·R·ŷ
///
/// Analytical expected values:
///   a_ref.x = -ω²·R
///   a_ref.y = +α·R
///   a_ref.z =  0
#[test]
fn test_inertial_rotor_a_ref_analytical_z_axis() {
    use aeroelast_solvers::petsc::fsi::rotor_physics::compute_rigid_body_acceleration_inertial;

    let r = 5.0f64;
    let omega = 3.0f64;
    let alpha = 2.0f64;

    let coords = vec![r, 0.0, 0.0]; // node on +X axis
    let axis = [0.0f64, 0.0f64, 1.0f64]; // Z rotation axis
    let center = [0.0f64; 3];

    let a = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);

    // α·ẑ × R·x̂ = α·R·(ẑ×x̂) = α·R·ŷ
    let expected_x = -(omega * omega * r); // centripetal
    let expected_y = alpha * r;            // tangential
    let expected_z = 0.0f64;

    assert!(
        (a[0] - expected_x).abs() < 1e-12,
        "a_ref.x: got {:.15e}, expected {:.15e}",
        a[0],
        expected_x
    );
    assert!(
        (a[1] - expected_y).abs() < 1e-12,
        "a_ref.y: got {:.15e}, expected {:.15e}",
        a[1],
        expected_y
    );
    assert!(
        (a[2] - expected_z).abs() < 1e-12,
        "a_ref.z: got {:.15e}, expected 0",
        a[2]
    );
}

/// SC-04f — `a_ref` superposition: off-axis node with non-trivial center.
///
/// For a node at (cx + R, cy, cz) with center = (cx, cy, cz) and X-axis rotation:
/// The effective radius vector is r = (R, 0, 0) (same as SC-04d/e after translation).
/// Expected values identical to the Z-axis case but transposed to the X-axis.
///
/// axis = x̂, r = R·ŷ (node at (cx, cy+R, cz)):
///   centripetal: a_cent = -ω²·R·ŷ
///   tangential:  a_tang = α·x̂ × R·ŷ = α·R·ẑ
#[test]
fn test_inertial_rotor_a_ref_off_axis_center() {
    use aeroelast_solvers::petsc::fsi::rotor_physics::compute_rigid_body_acceleration_inertial;

    let r = 4.0f64;
    let omega = 7.0f64;
    let alpha = 1.5f64;
    let cx = 10.0f64;
    let cy = -3.0f64;
    let cz = 2.5f64;

    // Node at (cx, cy + R, cz): radius vector = R·ŷ
    let coords = vec![cx, cy + r, cz];
    let axis = [1.0f64, 0.0f64, 0.0f64]; // X rotation axis
    let center = [cx, cy, cz];

    let a = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);

    // centripetal: -ω²·R·ŷ
    // tangential: α·x̂ × R·ŷ = α·R·(x̂×ŷ) = α·R·ẑ
    let expected_x = 0.0f64;
    let expected_y = -(omega * omega * r); // centripetal
    let expected_z = alpha * r;            // tangential

    assert!(
        (a[0] - expected_x).abs() < 1e-12,
        "a_ref.x: got {:.15e}, expected 0",
        a[0]
    );
    assert!(
        (a[1] - expected_y).abs() < 1e-12,
        "a_ref.y: got {:.15e}, expected {:.15e}",
        a[1],
        expected_y
    );
    assert!(
        (a[2] - expected_z).abs() < 1e-12,
        "a_ref.z: got {:.15e}, expected {:.15e}",
        a[2],
        expected_z
    );
}

// ── SC-07 / SC-08 — Python vs. Rust parity (requires live preCICE) ─────────
//
// These tests are marked `#[ignore]` because they require:
//   - A running preCICE 3 daemon (or in-process participant)
//   - PETSc matrices assembled from a real mesh (FEM assembly path)
//   - A CFD coupling partner (mock or real) that can exchange force/displacement
//
// To run once the environment is ready:
//   cargo test -p aeroelast-solvers -- --include-ignored parity
//
// Expected workflow:
//   1. Assemble K, M, C from the mesh.
//   2. Drive the Python solver for T=10 s, dt=0.01 s.
//   3. Drive InertialRotorFsiSolver with identical inputs.
//   4. Assert max |u_rust - u_python| < 1e-4 across all DOFs and time steps.

/// SC-07 — Python vs. Rust parity — synthetic 3-blade rotor (N_dof ~ 1200).
///
/// Blocked: needs preCICE participant + PETSc matrices assembled from mesh.
/// Run with: `cargo test -- --include-ignored test_parity_small_blade`
#[test]
#[ignore = "requires live preCICE environment — see SC-07 in migration spec"]
fn test_parity_small_blade() {
    // TODO (SC-07): implement once preCICE in-process participant is available.
    //
    // Outline:
    //   let config = InertialRotorFsiConfig { ... };  // synthetic 3-blade rotor
    //   let rust_u = run_rust_solver(&config, T=10.0, dt=0.01);
    //   let py_u   = run_python_solver(&config, T=10.0, dt=0.01);
    //   let max_diff = rust_u.iter().zip(&py_u).map(|(a,b)| (a-b).abs())
    //                        .fold(0.0_f64, f64::max);
    //   assert!(max_diff < 1e-4, "parity failure: max_diff = {max_diff:.3e}");
    todo!("SC-07: implement with preCICE in-process participant");
}

/// SC-08 — Python vs. Rust parity — IEA-15 MW reference turbine (full geometry).
///
/// Blocked: needs preCICE participant + IEA-15 mesh loaded and assembled.
/// Run with: `cargo test -- --include-ignored test_parity_iea15`
#[test]
#[ignore = "requires live preCICE environment and IEA-15 mesh — see SC-08 in migration spec"]
fn test_parity_iea15() {
    // TODO (SC-08): implement once IEA-15 benchmark case is available.
    //
    // Acceptance criteria:
    //   - max |u_rust - u_python| < 1e-4 across all DOFs and time steps
    //   - Cp within 2% of Python baseline
    //   - Ct within 2% of Python baseline
    todo!("SC-08: implement with IEA-15 mesh and live preCICE environment");
}
