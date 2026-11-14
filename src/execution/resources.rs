//! Admission accounting for total single-node budgets.
//!
//! Reservations are estimates used to decide what may start. They are not
//! operating-system limits; observed usage is recorded separately.

use serde::Serialize;

use crate::error::Error;

/// Total budgets shared by every concurrently running task.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Budget {
    pub threads: usize,
    pub memory_mib: u64,
    pub scratch_mib: u64,
    pub jobs: usize,
}

/// A helper process started alongside a tool, such as a sorter or compressor.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct Helper {
    pub name: String,
    pub threads: usize,
    pub memory_mib: u64,
}

/// What one task holds while it runs. Totals include every helper.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Reservation {
    pub threads: usize,
    pub memory_mib: u64,
    pub scratch_mib: u64,
    pub helpers: Vec<Helper>,
}

impl Reservation {
    pub fn total_threads(&self) -> usize {
        self.helpers.iter().fold(self.threads, |sum, helper| {
            sum.saturating_add(helper.threads)
        })
    }

    pub fn total_memory_mib(&self) -> u64 {
        self.helpers.iter().fold(self.memory_mib, |sum, helper| {
            sum.saturating_add(helper.memory_mib)
        })
    }
}

/// Resources currently reserved by running tasks.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Usage {
    pub threads: usize,
    pub memory_mib: u64,
    pub scratch_mib: u64,
    pub jobs: usize,
}

impl Budget {
    pub fn validate(&self) -> Result<(), Error> {
        if self.threads == 0 || self.jobs == 0 || self.memory_mib == 0 {
            return Err(Error::InvalidInput(
                "execution budget needs at least one thread, one job, and 1 MiB of memory".into(),
            ));
        }
        if self.jobs > self.threads {
            return Err(Error::InvalidInput(
                "concurrent jobs must not exceed the total thread budget".into(),
            ));
        }
        Ok(())
    }

    /// Rejects a task that could never be admitted, before anything launches.
    pub fn check(&self, task: &str, reservation: &Reservation) -> Result<(), Error> {
        let mut problems = Vec::new();
        if reservation.threads == 0 {
            problems.push("its main process reserves no thread".to_owned());
        }
        if reservation.total_threads() > self.threads {
            problems.push(format!(
                "{} threads including helpers exceed the {}-thread budget",
                reservation.total_threads(),
                self.threads
            ));
        }
        if reservation.total_memory_mib() > self.memory_mib {
            problems.push(format!(
                "{} MiB including helpers exceeds the {} MiB memory budget",
                reservation.total_memory_mib(),
                self.memory_mib
            ));
        }
        if reservation.scratch_mib > self.scratch_mib {
            problems.push(format!(
                "{} MiB of scratch exceeds the {} MiB scratch budget",
                reservation.scratch_mib, self.scratch_mib
            ));
        }
        if problems.is_empty() {
            Ok(())
        } else {
            Err(Error::InvalidInput(format!(
                "task {task} can never be admitted: {}",
                problems.join("; ")
            )))
        }
    }
}

/// Tracks reservations of running tasks against a budget.
#[derive(Debug)]
pub struct Ledger {
    budget: Budget,
    used: Usage,
    peak: Usage,
}

impl Ledger {
    pub fn new(budget: Budget) -> Self {
        Self {
            budget,
            used: Usage::default(),
            peak: Usage::default(),
        }
    }

    pub fn fits(&self, reservation: &Reservation) -> bool {
        self.used.jobs < self.budget.jobs
            && self
                .used
                .threads
                .saturating_add(reservation.total_threads())
                <= self.budget.threads
            && self
                .used
                .memory_mib
                .saturating_add(reservation.total_memory_mib())
                <= self.budget.memory_mib
            && self
                .used
                .scratch_mib
                .saturating_add(reservation.scratch_mib)
                <= self.budget.scratch_mib
    }

    /// Reserves resources if they fit; returns false and changes nothing otherwise.
    pub fn admit(&mut self, reservation: &Reservation) -> bool {
        if !self.fits(reservation) {
            return false;
        }
        self.used.jobs += 1;
        self.used.threads += reservation.total_threads();
        self.used.memory_mib += reservation.total_memory_mib();
        self.used.scratch_mib += reservation.scratch_mib;
        self.peak = Usage {
            threads: self.peak.threads.max(self.used.threads),
            memory_mib: self.peak.memory_mib.max(self.used.memory_mib),
            scratch_mib: self.peak.scratch_mib.max(self.used.scratch_mib),
            jobs: self.peak.jobs.max(self.used.jobs),
        };
        true
    }

    pub fn release(&mut self, reservation: &Reservation) {
        self.used.jobs -= 1;
        self.used.threads -= reservation.total_threads();
        self.used.memory_mib -= reservation.total_memory_mib();
        self.used.scratch_mib -= reservation.scratch_mib;
    }

    pub fn used(&self) -> Usage {
        self.used
    }

    pub fn peak(&self) -> Usage {
        self.peak
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn budget() -> Budget {
        Budget {
            threads: 4,
            memory_mib: 1000,
            scratch_mib: 100,
            jobs: 3,
        }
    }

    fn task(threads: usize, memory_mib: u64, helper_threads: usize) -> Reservation {
        Reservation {
            threads,
            memory_mib,
            scratch_mib: 10,
            helpers: vec![Helper {
                name: "sort".into(),
                threads: helper_threads,
                memory_mib: 100,
            }],
        }
    }

    #[test]
    fn helpers_count_toward_every_total() {
        let reservation = task(1, 200, 2);
        assert_eq!(reservation.total_threads(), 3);
        assert_eq!(reservation.total_memory_mib(), 300);
        let mut ledger = Ledger::new(budget());
        assert!(ledger.admit(&reservation));
        // One main thread is free, but the helper makes the second task need three.
        assert!(!ledger.admit(&reservation));
        assert_eq!(ledger.used().threads, 3);
    }

    #[test]
    fn admission_never_exceeds_any_budget_dimension() {
        let mut ledger = Ledger::new(budget());
        let small = task(1, 100, 0);
        for _ in 0..3 {
            assert!(ledger.admit(&small));
        }
        // Threads remain, but the job limit is reached.
        assert!(!ledger.fits(&small));
        ledger.release(&small);
        let memory_heavy = task(1, 900, 0);
        assert!(!ledger.admit(&memory_heavy));
        let scratch_heavy = Reservation {
            scratch_mib: 81,
            ..task(1, 1, 0)
        };
        assert!(!ledger.admit(&scratch_heavy));
        let used = ledger.used();
        assert!(used.threads <= 4 && used.memory_mib <= 1000 && used.scratch_mib <= 100);
        assert_eq!(ledger.peak().jobs, 3);
    }

    #[test]
    fn release_restores_capacity_exactly() {
        let mut ledger = Ledger::new(budget());
        let reservation = task(2, 500, 2);
        assert!(ledger.admit(&reservation));
        ledger.release(&reservation);
        assert_eq!(ledger.used(), Usage::default());
        assert!(ledger.admit(&reservation));
    }

    #[test]
    fn impossible_or_threadless_tasks_are_rejected_up_front() {
        let budget = budget();
        assert!(budget.check("ok", &task(1, 100, 1)).is_ok());
        for reservation in [
            task(0, 1, 0),
            task(4, 1, 1),
            task(1, 950, 0),
            Reservation {
                scratch_mib: 101,
                ..task(1, 1, 0)
            },
        ] {
            let error = budget.check("big", &reservation).unwrap_err();
            assert_eq!(error.exit_code(), 2);
            assert!(error.to_string().contains("task big can never be admitted"));
        }
        let overflow = Reservation {
            threads: usize::MAX,
            memory_mib: 1,
            scratch_mib: 0,
            helpers: vec![Helper {
                name: "h".into(),
                threads: 5,
                memory_mib: u64::MAX,
            }],
        };
        assert!(budget.check("overflow", &overflow).is_err());
    }

    #[test]
    fn budgets_need_capacity_and_jobs_within_threads() {
        assert!(budget().validate().is_ok());
        for invalid in [
            Budget {
                threads: 0,
                ..budget()
            },
            Budget {
                jobs: 0,
                ..budget()
            },
            Budget {
                memory_mib: 0,
                ..budget()
            },
            Budget {
                jobs: 5,
                ..budget()
            },
        ] {
            assert!(invalid.validate().is_err());
        }
    }
}
