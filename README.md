# Radial Dimensionality Theory IV: Tidal Deformability Ceilings under Solar-Constrained Density-Dependent Dimensionality

Scientific source and reproducible results for version v7.

- [Paper](reproduction/output/pdf/publication/paper4/merrill_rdt_tidal_deformability_v7.pdf)
- [Manuscript](reproduction/paper4_TidalLoveNumbers/V6_Paper/paper4_v6.tex)

## Reproduction

Requires Python 3.14, NumPy 2.4.2, SciPy 1.17.1, Matplotlib 3.10.8,
TeX Live with revtex4-2, and Poppler. The resource controller uses POSIX
nice/setpriority, SIGPROF and SIGSTOP/SIGCONT.

From the repository directory:

```sh
nice -n 19 python3 reproduce.py --output /tmp/merrill2025_tidal_deformability_code-reproduction
```

The output directory must not exist. This verifies every capsule checksum,
regenerates tables, figures or algebraic checks from retained scientific outputs,
and builds a PDF whose extracted text matches the included paper.
No parameter scan runs in this command. Full scientific runs use the producing
scripts and their data under `reproduction/`. Use one worker, verified nice 19,
at most 95% of one CPU and one native-library thread.

`release.json` identifies the current manuscript and PDF;
`reproduction/CAPSULE.json` binds the scientific source and data by checksum.
Shared dependencies and saved results are included so no sibling checkout is
required. Historical top-level scripts and data retain their original identity;
the current paper and source manifest distinguish them from the revised results.

## Scientific scope

This version retains the published title, but the word Ceilings no longer describes an established result: the revision reports conditional tidal samples, not global upper bounds. It integrates the full Love-number and radius response using the specified stellar backgrounds, GR-form perturbation equation and Schwarzschild exterior extraction. That exterior is a prescribed part of the calculation, not a matching result derived from the interior; the metric-gradient mismatch at the finite cutoff is about 1.30 percent in the primary 1.4-solar-mass example. This mismatch is a matching diagnostic, not an error bar on the tidal deformability.
