use std::path::Path;
use std::process::Command;

fn main() {
    // ── MPI ──────────────────────────────────────────────────────────────────
    // libpetsc.so depends on libmpi.so, which is usually pulled transitively by
    // PETSc's own pkg-config. We only emit explicit MPI flags when a concrete
    // install is found, so the build works in a conda env and on the HPC
    // cluster without hardcoding any absolute path.
    //
    // Resolution order:
    //   (a) conda env: <CONDA_PREFIX>/lib/libmpi.so
    //   (b) pkg-config for mpi / ompi / mpich
    //   (c) nothing (do not fail; PETSc's pkg-config may provide MPI transitively)
    let mpi_emitted = if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        let libmpi = Path::new(&prefix).join("lib").join("libmpi.so");
        if libmpi.exists() {
            let dir = libmpi.parent().unwrap().to_string_lossy().into_owned();
            println!("cargo:rustc-link-search=native={dir}");
            println!("cargo:rustc-link-lib=mpi");
            println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
            true
        } else {
            false
        }
    } else {
        false
    };

    if !mpi_emitted {
        for pkg in ["mpi", "ompi", "mpich"] {
            if let Some(flags) = pkg_config(pkg, &["--libs"]) {
                emit_link_flags(flags.as_bytes());
                emit_lib_flags(flags.as_bytes());
                if let Some(libdir) = pkg_config(pkg, &["--variable=libdir"]) {
                    println!("cargo:rustc-link-arg=-Wl,-rpath,{libdir}");
                }
                break;
            }
        }
    }

    // ── PETSc ────────────────────────────────────────────────────────────────
    let petsc_hint = "PETSc was not found via pkg-config. Activate a conda \
                      environment that contains petsc (e.g. `conda create -n \
                      aeroelast-dev -c conda-forge petsc slepc`) and rebuild.";
    let petsc_libs = require_pkg_config("PETSc", &["--libs-only-L"], petsc_hint);
    let petsc_link = require_pkg_config("PETSc", &["--libs-only-l"], petsc_hint);
    let petsc_rpath = require_pkg_config("PETSc", &["--variable=libdir"], petsc_hint);

    emit_link_flags(petsc_libs.as_bytes());
    emit_lib_flags(petsc_link.as_bytes());

    // rpath so the binary finds libpetsc.so at runtime
    if !petsc_rpath.is_empty() {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{petsc_rpath}");
    }

    // ── SLEPc ────────────────────────────────────────────────────────────────
    let slepc_hint = "SLEPc was not found via pkg-config. Activate a conda \
                      environment that contains slepc (e.g. `conda create -n \
                      aeroelast-dev -c conda-forge petsc slepc`) and rebuild.";
    let slepc_libs = require_pkg_config("SLEPc", &["--libs-only-L"], slepc_hint);
    let slepc_link = require_pkg_config("SLEPc", &["--libs-only-l"], slepc_hint);
    let slepc_rpath = require_pkg_config("SLEPc", &["--variable=libdir"], slepc_hint);

    emit_link_flags(slepc_libs.as_bytes());
    emit_lib_flags(slepc_link.as_bytes());

    // rpath so the binary finds libslepc.so at runtime
    if !slepc_rpath.is_empty() {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{slepc_rpath}");
    }

    // ── preCICE (feature = fsi) ───────────────────────────────────────────────
    // Only link if the `fsi` feature is active. Resolution order:
    //   (a) PRECICE_LIB_DIR env var (explicit override)
    //   (b) pkg-config --libs precice
    //   (c) conda env: <CONDA_PREFIX>/lib/libprecice.so
    //   (d) none of the above -> warn; the linker will fail with context
    if std::env::var("CARGO_FEATURE_FSI").is_ok() {
        let mut precice_found = false;

        if let Ok(dir) = std::env::var("PRECICE_LIB_DIR") {
            let dir = dir.trim().to_string();
            if !dir.is_empty() && Path::new(&dir).join("libprecice.so").exists() {
                println!("cargo:rustc-link-search=native={dir}");
                println!("cargo:rustc-link-lib=dylib=precice");
                println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
                precice_found = true;
            } else {
                println!(
                    "cargo:warning=PRECICE_LIB_DIR={dir} does not contain libprecice.so; falling back to other detection methods."
                );
            }
        }

        if !precice_found {
            if let Some(flags) = pkg_config("precice", &["--libs"]) {
                emit_link_flags(flags.as_bytes());
                emit_lib_flags(flags.as_bytes());
                if let Some(libdir) = pkg_config("precice", &["--variable=libdir"]) {
                    println!("cargo:rustc-link-arg=-Wl,-rpath,{libdir}");
                }
                precice_found = true;
            }
        }

        if !precice_found {
            if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
                let libprecice = Path::new(&prefix).join("lib").join("libprecice.so");
                if libprecice.exists() {
                    let dir = libprecice.parent().unwrap().to_string_lossy().into_owned();
                    println!("cargo:rustc-link-search=native={dir}");
                    println!("cargo:rustc-link-lib=dylib=precice");
                    println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
                    precice_found = true;
                }
            }
        }

        if !precice_found {
            println!(
                "cargo:warning=preCICE was not found (PRECICE_LIB_DIR, pkg-config, and CONDA_PREFIX all failed). The `fsi` feature will fail to link. Install preCICE with `conda install -c conda-forge precice pyprecice` and make sure the environment is active, or set PRECICE_LIB_DIR to the directory containing libprecice.so."
            );
        }

        println!("cargo:rerun-if-env-changed=PRECICE_LIB_DIR");
    }

    // Re-run build.rs if pkg-config output or the conda prefix changes
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_PATH");
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
}

/// Run `pkg-config <pkg> <args...>` and return trimmed stdout on success.
///
/// When a conda environment is active, prepend its pkg-config directories to
/// `PKG_CONFIG_PATH` for the subprocess so conda-installed libraries (PETSc,
/// SLEPc, ...) are found without requiring the user to set the variable.
fn pkg_config(pkg: &str, args: &[&str]) -> Option<String> {
    let mut cmd = Command::new("pkg-config");
    if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        let mut dirs = format!("{prefix}/lib/pkgconfig:{prefix}/share/pkgconfig");
        if let Ok(existing) = std::env::var("PKG_CONFIG_PATH") {
            if !existing.is_empty() {
                dirs.push(':');
                dirs.push_str(&existing);
            }
        }
        cmd.env("PKG_CONFIG_PATH", dirs);
    }
    let output = cmd.arg(pkg).args(args).output().ok()?;
    if !output.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&output.stdout).trim().to_string())
}

/// Like `pkg_config`, but aborts the build with a helpful message when the
/// package is missing or pkg-config itself cannot be run.
fn require_pkg_config(pkg: &str, args: &[&str], hint: &str) -> String {
    match pkg_config(pkg, args) {
        Some(out) => out,
        None => {
            println!("cargo:warning={hint}");
            panic!("pkg-config {pkg} {} failed. {hint}", args.join(" "));
        }
    }
}

/// Parse `-L/path/to/lib` flags and emit `cargo:rustc-link-search=native=...`
fn emit_link_flags(output: &[u8]) {
    let s = String::from_utf8_lossy(output);
    for token in s.split_whitespace() {
        if let Some(path) = token.strip_prefix("-L") {
            println!("cargo:rustc-link-search=native={path}");
        }
    }
}

/// Parse `-lfoo` flags and emit `cargo:rustc-link-lib=foo`
fn emit_lib_flags(output: &[u8]) {
    let s = String::from_utf8_lossy(output);
    for token in s.split_whitespace() {
        if let Some(lib) = token.strip_prefix("-l") {
            println!("cargo:rustc-link-lib={lib}");
        }
    }
}