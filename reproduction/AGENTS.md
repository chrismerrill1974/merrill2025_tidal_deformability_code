# Computational resource policy

The user requires one CPU per session so several AI sessions can run alongside each other.

- Run one computational worker only. Do not launch multiprocessing pools or
  concurrent scientific jobs.
- Establish and verify effective nice priority 19 before computation. A failed
  `nice`/`setpriority` call must stop the launch; do not assume it succeeded.
- Apply explicit CPU throttling to at most 95% of one core. Nice priority alone
  is not a CPU or thermal limit. Limit native numerical-library threads to one.
- Verify the live process priority, CPU use and absence of worker children.
- For step 15, use `paper3_NeutronStars/ns_step15/resource_control.py` as the
  launcher. Direct scientific entry points reject missing controls. Retain the
  resource-measurement logs with the run evidence.
- Preserve interrupted results and validate source/input compatibility and
  saved rows before reuse. Do not silently relabel a partial run complete.

These requirements implement the user's explicit resource constraints. They
apply to future research runs in this repository as well as step 15.
