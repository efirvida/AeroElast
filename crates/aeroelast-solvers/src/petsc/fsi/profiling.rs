/// Per-time-window profiling helper for FSI solvers.
///
/// Accumulates wall-clock timings of named code sections across the sub-iterations
/// of one FSI time window (preCICE implicit coupling can call the structural
/// solver several times per window before convergence). When the window
/// converges, [`WindowProfiler::log_summary`] emits a single `debug!` line with
/// the wall time of the whole window plus the **mean** cost of each tracked
/// section. The summary is silent when the `debug` log level is disabled, so
/// release builds carry only the cost of recording timings into a small `Vec`.
///
/// Usage pattern from a solver:
///
/// ```ignore
/// let mut prof = WindowProfiler::new("StressStiffened");
/// // start of converged window
/// prof.start_window();
/// // ... in the main loop, for each sub-iteration:
/// let t = Instant::now();
/// let forces = participant.read_data(...)?;
/// prof.record("precice_read", t.elapsed());
/// // ... and so on for every section of interest
/// // at the end of each converged window:
/// prof.log_summary(t_sim, window_idx);
/// ```
///
/// Sections are identified by `&'static str` for zero-alloc bookkeeping; an
/// internal `Vec` of `(name, total, count)` triplets is used instead of a
/// `HashMap` since the number of sections per solver is small (≲ 10) and
/// linear scan is faster at that size than hashing.
use std::time::{Duration, Instant};

/// A simple per-window profiler that buckets wall-clock timings by section name.
pub struct WindowProfiler {
    /// Display tag for log lines (typically the solver name).
    solver_tag: &'static str,
    /// Wall-clock instant when the current window started.
    window_start: Instant,
    /// Per-section (name, total accumulated time, number of samples) buckets.
    sections: Vec<(&'static str, Duration, u32)>,
}

impl WindowProfiler {
    /// Create a profiler tagged with the given solver name.
    pub fn new(solver_tag: &'static str) -> Self {
        Self {
            solver_tag,
            window_start: Instant::now(),
            sections: Vec::with_capacity(12),
        }
    }

    /// Reset accumulated samples and start the wall-clock timer for a new window.
    ///
    /// Call this once per converged FSI time window — typically right after the
    /// previous window's `log_summary`, or before the first sub-iteration of
    /// the first window.
    pub fn start_window(&mut self) {
        self.window_start = Instant::now();
        for entry in &mut self.sections {
            entry.1 = Duration::ZERO;
            entry.2 = 0;
        }
    }

    /// Record one timing sample for a named section.
    ///
    /// Multiple samples for the same name within a window are aggregated:
    /// totals are summed, counts incremented. The mean is computed at log time.
    pub fn record(&mut self, name: &'static str, elapsed: Duration) {
        if let Some(entry) = self.sections.iter_mut().find(|e| e.0 == name) {
            entry.1 += elapsed;
            entry.2 += 1;
        } else {
            self.sections.push((name, elapsed, 1));
        }
    }

    /// Convenience: time the closure `f` and record the elapsed duration under `name`.
    ///
    /// Useful when the timed block is a single expression and doesn't need to
    /// use `?` early-return. For blocks that propagate errors with `?`, use
    /// `record_since(name, t)` with an `Instant` taken manually.
    pub fn time<F, R>(&mut self, name: &'static str, f: F) -> R
    where
        F: FnOnce() -> R,
    {
        let t = Instant::now();
        let r = f();
        self.record(name, t.elapsed());
        r
    }

    /// Convenience: record the time elapsed since a previously-captured `Instant`.
    ///
    /// Intended for sections that use `?` for error propagation:
    ///
    /// ```ignore
    /// let t = Instant::now();
    /// participant.read_data(...)?;
    /// prof.record_since("precice_read", t);
    /// ```
    pub fn record_since(&mut self, name: &'static str, since: Instant) {
        self.record(name, since.elapsed());
    }

    /// Emit the per-window summary line at `debug` level.
    ///
    /// Format:
    /// `{solver_tag} window {idx} t={t_sim}s wall={wall}ms | {sec1}: {mean1}ms (×{n1}), ...`
    ///
    /// `mean_ms` is `total_ms / count` so a section called twice gets its
    /// per-call average reported. The wall time is the elapsed time since
    /// `start_window` and does NOT divide by anything — it is the total
    /// wall-clock cost of the window (preCICE coupling iterations included).
    ///
    /// Silent when the `debug` level is disabled; the only overhead in release
    /// is the level check and the prior `record` accumulations.
    pub fn log_summary(&self, t_sim: f64, window_idx: usize) {
        if !log::log_enabled!(log::Level::Debug) {
            return;
        }
        let wall_ms = self.window_start.elapsed().as_secs_f64() * 1000.0;
        let mut msg = format!(
            "{} window {} t={:.4}s wall={:.1}ms",
            self.solver_tag, window_idx, t_sim, wall_ms
        );
        for (name, total, count) in &self.sections {
            if *count == 0 {
                continue;
            }
            let mean_ms = total.as_secs_f64() * 1000.0 / (*count as f64);
            msg.push_str(&format!(", {}: {:.2}ms (×{})", name, mean_ms, count));
        }
        log::debug!("{}", msg);
    }

    /// Test/diagnostic accessor: returns the wall-clock cost of the current window so far.
    #[cfg(test)]
    pub fn current_wall(&self) -> Duration {
        self.window_start.elapsed()
    }

    /// Test/diagnostic accessor: total accumulated time for a section in the current window.
    #[cfg(test)]
    pub fn section_total(&self, name: &'static str) -> Option<Duration> {
        self.sections.iter().find(|e| e.0 == name).map(|e| e.1)
    }

    /// Test/diagnostic accessor: sample count for a section in the current window.
    #[cfg(test)]
    pub fn section_count(&self, name: &'static str) -> Option<u32> {
        self.sections.iter().find(|e| e.0 == name).map(|e| e.2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread::sleep;

    #[test]
    fn test_record_accumulates_across_samples() {
        let mut p = WindowProfiler::new("Test");
        p.start_window();
        p.record("foo", Duration::from_millis(10));
        p.record("foo", Duration::from_millis(20));
        p.record("bar", Duration::from_millis(5));

        assert_eq!(p.section_count("foo"), Some(2));
        assert_eq!(p.section_count("bar"), Some(1));
        assert_eq!(p.section_total("foo"), Some(Duration::from_millis(30)));
    }

    #[test]
    fn test_start_window_clears_samples() {
        let mut p = WindowProfiler::new("Test");
        p.start_window();
        p.record("foo", Duration::from_millis(10));
        assert_eq!(p.section_count("foo"), Some(1));

        p.start_window();
        // Sample count reset, but the section slot is retained (allocation reuse).
        assert_eq!(p.section_count("foo"), Some(0));
        assert_eq!(p.section_total("foo"), Some(Duration::ZERO));
    }

    #[test]
    fn test_time_closure_records_elapsed() {
        let mut p = WindowProfiler::new("Test");
        p.start_window();
        let result = p.time("compute", || {
            sleep(Duration::from_millis(5));
            42
        });
        assert_eq!(result, 42);
        assert_eq!(p.section_count("compute"), Some(1));
        // Recorded time is at least the sleep duration.
        assert!(p.section_total("compute").unwrap() >= Duration::from_millis(5));
    }

    #[test]
    fn test_record_since_uses_instant() {
        let mut p = WindowProfiler::new("Test");
        p.start_window();
        let t = Instant::now();
        sleep(Duration::from_millis(3));
        p.record_since("io", t);
        assert_eq!(p.section_count("io"), Some(1));
        assert!(p.section_total("io").unwrap() >= Duration::from_millis(3));
    }

    #[test]
    fn test_window_wall_clock_grows() {
        let mut p = WindowProfiler::new("Test");
        p.start_window();
        let w0 = p.current_wall();
        sleep(Duration::from_millis(2));
        let w1 = p.current_wall();
        assert!(w1 > w0);
    }

    #[test]
    fn test_log_summary_silent_when_debug_disabled() {
        // Default test harness has no env_logger initialized → debug disabled.
        // This test mostly confirms log_summary is a no-panic operation in that case.
        let mut p = WindowProfiler::new("Test");
        p.start_window();
        p.record("x", Duration::from_millis(1));
        p.log_summary(1.234, 42); // must not panic regardless of log level
    }
}
