/// FSI (Fluid-Structure Interaction) coupling module.
///
/// Provides the building blocks for preCICE-based FSI coupling:
/// - `checkpoint` — state snapshot for implicit coupling iteration rollback
/// - `force_utils` — force ramping and capping before passing to the stepper
/// - `profiling`   — `WindowProfiler` for per-time-window component timing breakdown
/// - `setup`       — BC reduction, mass lumping, Rayleigh C and interface DOF remapping
/// - `linear_elastic` — full `LinearElasticFsiSolver` coupling loop
/// - `stress_stiffened` — `StressStiffenedFsiSolver` with per-step K_G update
/// - `rotor_physics` — coordinate transforms, inertial forces, K_SP, torque, OmegaProvider
/// - `rotor_fsi`     — `RotorFsiSolver` full co-rotational FSI coupling loop
/// - `rotor_inertial` — `InertialRotorFsiSolver` full inertial-frame FSI coupling loop

pub mod checkpoint;
pub mod force_utils;
pub mod linear_elastic;
pub mod profiling;
pub mod rotor_fsi;
pub mod rotor_inertial;
pub mod rotor_physics;
pub mod setup;
pub mod stress_stiffened;

pub use linear_elastic::{FsiConfig, FsiError, FsiInitialState, FsiResult, LinearElasticFsiSolver};
pub use rotor_fsi::RotorFsiSolver;
pub use rotor_inertial::InertialRotorFsiSolver;
pub use stress_stiffened::StressStiffenedFsiSolver;
