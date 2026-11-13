# One CPU per session

The user permits this session up to 95% of one CPU, at nice priority 19, so
other AI sessions can run their own jobs. One computational worker and one
native numerical thread are required. Do not run scientific stages concurrently.

`resource_control.py` establishes and checks priority, sets six numerical-library
thread environment variables to one before numerical imports, and enforces a
CPU-time budget using `ITIMER_PROF`. Each 10 ms process-CPU quantum is followed
by any required sleep; idle time does not accumulate spendable CPU credit.
The scientific entry points check that the launcher is active. This is a soft
measured process throttle; scheduler and signal-delivery granularity can affect
short observation windows. Nice priority alone is not the CPU cap.

`resource_evidence/probe95.json` records a successful 2-second CPU probe:
2.000219 CPU seconds over 2.153851 wall seconds (92.87%), nice 19 and 188
throttle sleeps. `live95_process.txt` records the resumed process at 93.3% CPU,
nice 19; `child95_pids.txt` is empty. Earlier 25% evidence is preserved separately.
Each stage's resource log records final CPU time, wall time, limits and status.
The initial running-log measurement spans microseconds and is not a meaningful
utilization estimate; use the probe, live observation and final measurement.

The host process session allows short monitoring calls to return while the
calculation continues. The scan flushes a complete JSON journal row after each
parameter point. There is no overall scientific-run deadline configured.
Monitoring is available while this assistant session remains active; no promise
is made of monitoring after application/session termination. Throttling is
implemented within the process and does not depend on monitoring calls.

Never launch archived parallel scripts. They are retained only as provenance.
