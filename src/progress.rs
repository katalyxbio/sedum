use std::io::{stderr, IsTerminal, Write};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

/// Shared atomic counters updated by scan threads and read by a reporter
/// thread. All reads/writes are `Relaxed` — progress display does not need
/// any ordering guarantees.
pub struct Progress {
    pub stripes_total: AtomicU64,
    pub stripes_done: AtomicU64,
    pub records: AtomicU64,
    stop: AtomicBool,
    start: Instant,
    enabled: bool,
}

impl Progress {
    pub fn new(enabled: bool) -> Arc<Self> {
        Arc::new(Self {
            stripes_total: AtomicU64::new(0),
            stripes_done: AtomicU64::new(0),
            records: AtomicU64::new(0),
            stop: AtomicBool::new(false),
            start: Instant::now(),
            enabled,
        })
    }

    pub fn set_stripes_total(&self, n: u64) {
        self.stripes_total.store(n, Ordering::Relaxed);
    }

    #[inline]
    pub fn inc_stripe(&self) {
        self.stripes_done.fetch_add(1, Ordering::Relaxed);
    }

    #[inline]
    pub fn inc_records(&self, n: u64) {
        if n > 0 {
            self.records.fetch_add(n, Ordering::Relaxed);
        }
    }
}

/// RAII guard: spawns a reporter thread on `new`, joins and prints the
/// final summary on `drop`.
pub struct Reporter {
    handle: Option<JoinHandle<()>>,
    progress: Arc<Progress>,
}

impl Reporter {
    pub fn spawn(progress: Arc<Progress>) -> Self {
        let handle = if progress.enabled {
            let p = Arc::clone(&progress);
            Some(thread::spawn(move || run(p)))
        } else {
            None
        };
        Self { handle, progress }
    }
}

impl Drop for Reporter {
    fn drop(&mut self) {
        self.progress.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            // Wake the reporter immediately instead of waiting out its current
            // sleep interval (up to 5s when stderr is not a TTY).
            h.thread().unpark();
            let _ = h.join();
        }
        if !self.progress.enabled {
            return;
        }
        let p = &self.progress;
        let elapsed = p.start.elapsed().as_secs_f64().max(1e-9);
        let records = p.records.load(Ordering::Relaxed);
        let stripes_done = p.stripes_done.load(Ordering::Relaxed);
        let stripes_total = p.stripes_total.load(Ordering::Relaxed);
        let rate = records as f64 / elapsed;
        let mut err = stderr().lock();
        if stripes_total > 0 {
            let _ = writeln!(
                err,
                "\rsedum: done {stripes_done}/{stripes_total} stripes, \
                 {records} records in {elapsed:.1}s ({rate:.0} rec/s)          "
            );
        } else {
            let _ = writeln!(
                err,
                "\rsedum: done {records} records in {elapsed:.1}s ({rate:.0} rec/s)          "
            );
        }
    }
}

fn run(p: Arc<Progress>) {
    let tty = stderr().is_terminal();
    let interval = if tty {
        Duration::from_millis(250)
    } else {
        Duration::from_secs(5)
    };
    let cr = if tty { '\r' } else { '\n' };
    // Brief initial delay so we don't print before any work has happened
    // (park_timeout so a quick run can wake us out of it on shutdown).
    thread::park_timeout(Duration::from_millis(100));
    while !p.stop.load(Ordering::Relaxed) {
        let elapsed = p.start.elapsed().as_secs_f64().max(1e-9);
        let records = p.records.load(Ordering::Relaxed);
        let rate = records as f64 / elapsed;
        let stripes_total = p.stripes_total.load(Ordering::Relaxed);
        let stripes_done = p.stripes_done.load(Ordering::Relaxed);
        let mut err = stderr().lock();
        if stripes_total > 0 {
            let pct = 100.0 * stripes_done as f64 / stripes_total as f64;
            let eta = if stripes_done > 0 {
                let per_stripe = elapsed / stripes_done as f64;
                let remaining = (stripes_total - stripes_done) as f64 * per_stripe;
                format!("{remaining:>5.0}s")
            } else {
                "  --s".to_string()
            };
            let _ = write!(
                err,
                "{cr}sedum: {pct:5.1}% [{stripes_done}/{stripes_total} stripes] \
                 {records:>12} records  {rate:>9.0} rec/s  eta {eta}    "
            );
        } else {
            let _ = write!(
                err,
                "{cr}sedum: {records:>12} records  {rate:>9.0} rec/s  elapsed {elapsed:>6.1}s    "
            );
        }
        let _ = err.flush();
        drop(err);
        // park_timeout instead of sleep so `drop` can wake us instantly via
        // unpark; a spurious wakeup just re-checks `stop` and reprints.
        thread::park_timeout(interval);
    }
}
