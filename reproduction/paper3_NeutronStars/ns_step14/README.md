# Step 14 neutron-star numerical audit

Versioned working audit; historical solvers/results are preserved.
See [the report](../../2026_june_14_make_it_consistent/STEP14_NS_NUMERICS.md)
for equations, results and limitations. Computational collaborator: Astra (OpenAI).

- `model.py`: independent pressure-coordinate flux solver, checked EOS loader,
  corrected central series; linear-log or log-PCHIP EOS interpolation.
- `direct.py`: independent radial-field solver, exact interpolant derivative
  or legacy-stencil diagnostic; checked terminal surface.
- `run.py`: all 11 paired settings, source-convention regression, selected
  independent comparisons; every evaluation journaled with a completion manifest.
- `verify.py`: verify hashes, bind summary outputs to journal coordinates,
  recompute fractional shifts and measured sensitivity ranges.
- `test_model.py`: failure-path, central-series, GR analytic and independent
  implementation tests (11 methods, with parameterized subtests).
- `fixed_mass_crosscheck.py`: radial-field fixed-mass roots, with recorded
  settings; default and tightened output files retain both comparisons.
- `audit_sources.py`: compare official ZIP members and original DH PDF table
  values with repository inputs. Requires Poppler `pdftotext`.
- `endpoint_probe.py`: resolve the shallow GR turnover and distinguish it from
  the increasing modified curves at the available EOS edge.
- `reproduce.py`: stage the complete numerical source/input set in a fresh
  directory and repeat tests, calculations, verification and fixed-mass checks;
  compare results byte-for-byte with this revision.

From repository root (destination must not exist):

```sh
PYTHONDONTWRITEBYTECODE=1 python3 paper3_NeutronStars/ns_step14/reproduce.py --output-root /tmp/new-rdt-step14
```

Python dependency versions are pinned in `requirements.txt`; the recorded run
used Python 3.14.3. The driver uses the current interpreter. It does not install
packages, download tables, or silently reuse previous run outputs.

To repeat source verification, download the files explicitly:

```sh
curl -fL https://compose.obspm.fr/download/1D/NS/Skyrme/SLY4/eos.zip -o /tmp/compose-sly4.zip
curl -fL https://arxiv.org/pdf/astro-ph/0111092 -o /tmp/dh2001.pdf
PYTHONDONTWRITEBYTECODE=1 python3 paper3_NeutronStars/ns_step14/audit_sources.py --compose-zip /tmp/compose-sly4.zip --dh-pdf /tmp/dh2001.pdf --output /tmp/rdt-source-check.json
```

Recorded URLs/checksums and the detailed match results are in `source_audit.json`.
A future upstream revision is not automatically adopted into this study.

Results are at the historical test point A=.0087, alpha=.12, with a few explicit
cross-check probes; they are not the step 15 solar-domain search. The modified
RG(SLy4) masses reach the EOS boundary and are not established maximum masses;
GR has a shallow interior turnover resolved by near-edge refinement.
The surface is the finite tabulated pressure floor, not an extrapolated vacuum.
Measured sensitivity ranges are not statistical or rigorous error bounds.

The superseded 18-point initial grid and its source snapshots are retained in
`diagnostics/initial_grid`; current results use 22 grid points plus refinement.
Historical source snapshots are for restoration to their original paths in a
fresh checkout, not execution from the archive directory.
