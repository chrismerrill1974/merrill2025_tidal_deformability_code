# Solar parameter-region recomputation (step 11)

Computational collaborator: Astra (OpenAI).

This is a separate successor calculation. `solarcore_v3/results_scan.csv`, its
figures and the released papers are historical inputs, not overwritten outputs.
The scientific specification is `STEP09_INFERENCE_CONTRACT.json`; the retained
physical response is `STEP08_RESPONSE_DECISION.json` in the review directory.

`solver.py` integrates the same regularized expanded n=3 operator with ordinary
three-dimensional mass and reaction kernels as additional ODE states. Their
central integrals start at xi0^3/3 (neglected corrections are O(xi0^5)). The
calibration still requires relative density and reference-log residuals below
1e-12, with at most 60 iterations. Profiles use DOP853 dense output rather than
linear interpolation on an adaptive mesh. No EOS/composition/luminosity response
or measured helioseismic likelihood has been added.

`inference.py` implements the source-precision score and all named sensitivity
variants. Every fixed-alpha set uses its own minimum and DeltaQ=2.71. All sign
changes on the checked mesh are refined using actual model evaluations;
connected accepted intervals, including grid-limited intervals, are retained.
The code explicitly checks for sub-grid minima next to the A=0 boundary. Mesh
checks are numerical evidence, not a proof that arbitrarily narrow unsampled
components cannot exist. Alpha envelopes are sampled, not certified global
bounds or simultaneous confidence bands.

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_step11_inference.py
PYTHONDONTWRITEBYTECODE=1 python3 paper1_SolarCore/solar_step11/produce.py --output /tmp/new-solar-a
PYTHONDONTWRITEBYTECODE=1 python3 paper1_SolarCore/solar_step11/reference_run.py /tmp/new-solar-R
```

`produce.py` is the supported entry point: it requires the numerical/inference
checks to pass before invoking the runner. The recorded runs were launched
with the same tests as an explicit preceding command. Output directories must
be new. There is no automatic cache adoption or resume. A run failure leaves
its manifest failed and its evaluated points visible. Analysis requires a
complete manifest, matching source/input hashes, checksums, declared grid,
unique finite rows and converged calibration. These are single-writer artifacts.

`reference_run.py` explicitly selects r_star=R_sun and the tighter ODE settings
(rtol=2e-13, atol=2e-15, maximum step .02), after the preserved initial reference
run failed calibration. It also runs the gates first.

`audit.py` performs the separately recorded solver, units, reference-start,
profile-grid and surface-floor comparisons. `local_checks.py` preserves failed
wide-bracket checks and retries near known baseline endpoints, recording the
amplitude of any new failure. `summarize.py` admits only complete primary/reference
runs and calculates comparison summaries. No released figure is regenerated here.

The derivative cutoff remains a numerical prescription. Failed stricter settings
are retained; successful comparisons do not certify the zero-cutoff limit. The
legacy scorer uses exactly the earlier rounded inputs and selective theory
widths, evaluated on the new response; it is not a replacement for archived
old response values. No statistic here establishes a detection, confidence
coverage, or a realistic solar-model uncertainty.

## Reproduce in a clean directory

Install the recorded dependencies from `requirements.txt`; the checked runtime
was Python 3.14.3 (the complete runtime string is in the manifests). Then:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 paper1_SolarCore/solar_step11/reproduce.py --output-root /tmp/new-rdt-study
```

The destination must not exist. This copies required sources and input files,
records a source manifest, runs the gates and both scans, then runs the numerical
audits, local alpha refinements, summary and verification. It captures command
logs and stops on a failing command. This is a full numerical study, not an
instant report-generation command. Add `--stage-only` to copy/hash the sources
without running it. Historical failed development attempts are retained in this
repository for provenance; the reproduction driver executes the final workflow.

