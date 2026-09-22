/// RAII wrappers for opaque PETSc/SLEPc handles that are created, configured,
/// and destroyed within a single scope (KSP, EPS, SNES).
///
/// Each wrapper owns its handle and calls the matching `*Destroy` on drop, so
/// a handle is never leaked even if a configuration call fails mid-way.
use super::ffi::{self, EPS, KSP, SNES};
use std::fmt;

/// Safe wrapper around a PETSc `KSP` handle.
pub struct PetscKsp {
    pub(crate) raw: KSP,
}

impl PetscKsp {
    /// Wrap a raw PETSc KSP handle (takes ownership).
    ///
    /// # Safety
    /// `raw` must be a valid, initialized PETSc KSP that has not yet been destroyed.
    pub unsafe fn from_raw(raw: KSP) -> Self {
        Self { raw }
    }

    /// Return the raw handle (for passing to FFI calls).
    pub fn as_raw(&self) -> KSP {
        self.raw
    }
}

impl Drop for PetscKsp {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            // SAFETY: self.raw is a valid non-null KSP owned by this wrapper.
            unsafe {
                ffi::KSPDestroy(&mut self.raw);
            }
        }
    }
}

impl fmt::Debug for PetscKsp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PetscKsp({:p})", self.raw)
    }
}

// PETSc handles can be passed between threads as long as there is no
// concurrent access (caller's responsibility).
unsafe impl Send for PetscKsp {}

/// Safe wrapper around a SLEPc `EPS` handle.
pub struct PetscEps {
    pub(crate) raw: EPS,
}

impl PetscEps {
    /// Wrap a raw SLEPc EPS handle (takes ownership).
    ///
    /// # Safety
    /// `raw` must be a valid, initialized SLEPc EPS that has not yet been destroyed.
    pub unsafe fn from_raw(raw: EPS) -> Self {
        Self { raw }
    }

    /// Return the raw handle (for passing to FFI calls).
    pub fn as_raw(&self) -> EPS {
        self.raw
    }
}

impl Drop for PetscEps {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            // SAFETY: self.raw is a valid non-null EPS owned by this wrapper.
            unsafe {
                ffi::EPSDestroy(&mut self.raw);
            }
        }
    }
}

impl fmt::Debug for PetscEps {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PetscEps({:p})", self.raw)
    }
}

// PETSc handles can be passed between threads as long as there is no
// concurrent access (caller's responsibility).
unsafe impl Send for PetscEps {}

/// Safe wrapper around a PETSc `SNES` handle.
pub struct PetscSnes {
    pub(crate) raw: SNES,
}

impl PetscSnes {
    /// Wrap a raw PETSc SNES handle (takes ownership).
    ///
    /// # Safety
    /// `raw` must be a valid, initialized PETSc SNES that has not yet been destroyed.
    pub unsafe fn from_raw(raw: SNES) -> Self {
        Self { raw }
    }

    /// Return the raw handle (for passing to FFI calls).
    pub fn as_raw(&self) -> SNES {
        self.raw
    }
}

impl Drop for PetscSnes {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            // SAFETY: self.raw is a valid non-null SNES owned by this wrapper.
            unsafe {
                ffi::SNESDestroy(&mut self.raw);
            }
        }
    }
}

impl fmt::Debug for PetscSnes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PetscSnes({:p})", self.raw)
    }
}

// PETSc handles can be passed between threads as long as there is no
// concurrent access (caller's responsibility).
unsafe impl Send for PetscSnes {}