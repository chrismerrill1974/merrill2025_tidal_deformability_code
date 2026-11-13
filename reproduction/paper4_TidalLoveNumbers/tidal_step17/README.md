# Step 17 — conditional tidal-response audit

This versioned calculation retains the full density law, both checked background
prescriptions and the CompOSE RG(SLy4) EOS. It integrates the quadrupole response;
it does not infer tides from radius scaling. The GR-form perturbation equation,
Schwarzschild Love extraction and finite-cutoff matching are explicitly assumed.
See `2026_june_14_make_it_consistent/STEP17_TIDAL_RESPONSE.md` for equations,
controls, results and the measured exterior-gradient mismatch.

All stages run **one process/thread, nice 19, max 95% of one CPU**. Use `launch.py`;
direct computational functions refuse missing controls. Do not run stages in
parallel. The launcher sets limits before importing numerical libraries and
writes final resource measurements even on failure.

From the repository root, examples for fresh output directories:

```sh
nice -n 19 python3 paper4_TidalLoveNumbers/tidal_step17/launch.py \
  --log /tmp/step17-gates-resource.json \
  paper4_TidalLoveNumbers/tidal_step17/gates.py --output /tmp/step17-new-gates
nice -n 19 python3 paper4_TidalLoveNumbers/tidal_step17/launch.py \
  --log /tmp/step17-run-resource.json \
  paper4_TidalLoveNumbers/tidal_step17/run.py --output /tmp/step17-new-results
```

The run checks the committed admission-gate record and step 15 source/output
contract, then hashes its own sources and gates in its new manifest. Compare
fresh gates against `gates/gates.json`; do not replace failed or historical
artifacts. `verify.py` recomputes bindings and identities for the retained
`results` directory; `plot.py` plots those verified results into a fresh
`figures` directory. Their output paths are deliberately fixed to this retained
package. A full fresh repeat of the 60-star suite was not performed in this
step; selected paired tightening and alternative-formulation tests were.

- `tides.py`: EOS-knot-resolved pressure-coordinate solver, exact interpolant
  derivative, corrected central seeds, Riccati and linear-H forms, explicit
  matching alternatives and exact response decomposition.
- `gates.py`: Newtonian analytic checks, GR/background checks, energy-derivative
  check, invalid-input rejection and injected integration failure.
- `run.py`: 60 retained stellar solves, 21 sampled pairs, 24 paired controls.
- `verify.py`: parameter/output bindings, extraction algebra, pairing and binaries.
- `legacy_replay.py`: unchanged historical V6 example, separate from current output.
- `plot.py`: retained figure generator.
- `launch.py`: mandatory resource controls and measurements.

`results/evaluations.jsonl` stores every full stellar output, all pressure-knot
endpoints, all matching alternatives, and all root-search background evaluations.
Each completed stellar row is flushed immediately. Failure cannot mark a run
complete. Source and artifact hashes are in the run and step-level manifests.
