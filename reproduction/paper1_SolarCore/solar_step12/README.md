# Step 12 solar figures and table

Computational collaborator: Astra (OpenAI).

`figures.py` admits the complete step 11 solar-a and solar-R artifacts, verifies
source/input hashes, and generates score maps, slice endpoints and the selected
endpoint table. Contours are nominal score increments, not sigma levels.
Grid-limited endpoints remain unresolved. The summary also records finite-core
and asymptotic dimensions and exact inherited fiducial examples.

All generated files are in `figures/`, with a manifest. Original figures remain
historical. Paper I is built with these files by the retained
`2026_june_14_make_it_consistent/build_step12_papers.py` script, which stages
current figure/table inputs in a temporary directory.

To reproduce Paper I from this archive, run the repository-root `reproduce.py`.
It verifies the source/data hashes and regenerates the solar figures and tables.
Full solar recomputation uses the step 11 driver and the required CPU controls.
