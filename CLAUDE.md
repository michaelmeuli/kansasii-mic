# CLAUDE.md

See README.md for the question, data handling, point-estimate conventions and
outputs. This file only lists what is easy to get wrong.

- **Git:** commit and push directly to `main` (remote `origin`,
  `git@github.com:michaelmeuli/kansasii-mic.git`). No feature branches or PRs.
- **Environment:** the system `python3` has no pandas. Use
  `module load miniforge3` and `conda activate kansasii_mic`. Run tests with
  `pytest tests/`. Set `MPLBACKEND=Agg` when running headless.
- **Screening map is `data/imm/screening_map_link.csv`** (built by
  `immensekansasii/scripts/screening_map_link.py`); the old `screening_map.csv`
  is superseded. A `mic.csv` TNR matches any of `SCREENING_MAP_TNR_COLUMNS`
  (`TNR`, `TNR_NGS`, `TNR3`-`TNR6`), then the map's `MHK` column. The output
  `TNR` is always the isolate's primary `TNR`. The map's `MHK` column is
  unrelated to `mic.csv`'s `MHK` (the MIC reading).
- **Outputs are on the shares, not in the repo:**
  `/shares/sander.imm.uzh/MM/kansasii/output/mic/` (`mhk/` and `mgit/`).
  Before a rerun, empty that directory first so no stale figures remain.
  Inputs are in `/shares/sander.imm.uzh/MM/kansasii/data/imm/`.
- **Two parallel pipelines** share the same screening-map matching:
  MHK broth-microdilution rows (numeric MIC) and MGIT breakpoint rows
  (labels like `"Amikacin 4 mg/l"`, categorical, empty MHK). Don't mix them.
- **`scripts/run.sh`** is a Windows/PowerShell scratch log (uv venv, scp to the
  cluster), not a runnable bash script.
