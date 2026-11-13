# Step 12 white-dwarf propagation

Computational collaborator: Astra (OpenAI).

`propagate.py` reads the verified, pinned step 11 solar-a sets and local peak,
then calls the step 07 audited flux solver. It samples 66 boundary/exponent
points, five amplitude fractions and 12 densities. It requires a new output
directory, rejects incompatible solar input and records any failed run.

`figures.py` reads complete, hash-checked output and generates the actual
comparison-mass plot and exact fractional-shift table in a new directory.
The numerical outputs are samples, not a certified continuum envelope. The
upper density range exceeds the non-relativistic EOS applicability estimate.
The historical derivative-cutoff solver failures and archive discrepancy remain
in their original audit; this run explicitly uses the audited flux formulation.

Run the repository-root `reproduce.py` to regenerate the figure and table
from the retained results and build Paper II. Full white-dwarf propagation is
performed by `propagate.py`, using the pinned solar inputs and audited flux
solver. Apply the required one-worker, nice-19 and single-CPU throttle controls
to all scientific computation. Saved rows and source/input bindings are retained
in `results/manifest.json` and the associated numerical evidence.
