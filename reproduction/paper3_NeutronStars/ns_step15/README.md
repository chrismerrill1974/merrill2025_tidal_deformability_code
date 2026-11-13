# Step 15 — solar-region dependence of conditional NS prescriptions

All scientific work must use **one worker, verified nice 19, and explicit CPU
throttling**. The retained launcher limits process CPU to 95% of one core in
10-ms CPU quanta, disables accumulated idle credit and limits native numerical
threads to one. This is a measured process limit, not a guarantee about hardware
temperature or other applications. See the repository `AGENTS.md`.

Never launch the archived four-worker scripts. `diagnostics/scan_interrupted.py`
is retained solely to establish the source provenance of saved results.

Example launch from repository root:

```sh
nice -n 19 python3 paper3_NeutronStars/ns_step15/resource_control.py \
  --log /tmp/step15-resource.json \
  paper3_NeutronStars/ns_step15/scan.py --output /tmp/new-step15-scan --workers 1
```

If the environment refuses priority changes, resolve that before computation;
never continue at normal priority. The launcher itself checks priority and
throttling, and the numerical entry points refuse missing controls.

The original interrupted journal is under `results`. The completed run is under
`resumed95_results`, with an adoption record verifying source hashes, unchanged
numerical function AST, grid membership, result binding, numerical algebra and
exact baseline/selected-row replay. The original files remain unchanged.

The scan covers 66 solar parameter slices (65 saved slices and the locally
refined solar-amplitude peak), five amplitude fractions, both prescriptions,
and physical reference lengths of 5, 10 and 20 km. These are finite samples,
not a continuum confidence envelope. EOS-table and finite-surface limitations
from step 14 continue to apply.

Producing scripts:

- `resource_control.py`: mandatory launcher and measured CPU-budget probe.
- `scan.py`: fresh serial scan using the step 14 solver.
- `resume.py`: validate/adopt the interrupted journal and calculate missing rows.
- `refine.py`: recompute solar endpoints and refine observed NS response extrema.
- `checks.py`: numerical controls, independent field roots and saturation tests.
- `analyze.py`: bind outputs to input slices and recompute summaries.
- `plot.py`: generate the verified parameter-dependence figure.
- `reproduce.sh`: sequential fresh-directory reproduction (full repeat not run here).
- `verify_artifacts.py`: check the retained artifact manifest and completion evidence.
- `baseline_check.py`: retain the initial failed-precision gate and its resolution.

Historical diagnostic scripts are snapshots, not current execution commands.

The initial controlled run used 25% of one CPU; the continuation used a 95%
single-CPU limit. Its 1,597 completed rows are retained in `resumed_results`
and adopted into `resumed95_results`; archived 25% sources are in
`diagnostics/checkpoint25`. Other AI sessions have their own CPU budgets.
