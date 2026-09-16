# kansasii-mic

Outlier analysis of *Mycobacterium kansasii* antimicrobial MIC (minimum
inhibitory concentration) data: which isolates (`TNR`) are unusually
resistant or susceptible relative to the rest of the cohort, per antibiotic
and overall, with results checked against CLSI reference breakpoints for
*M. kansasii*.

## Data

- `data/imm/mic.csv` -- one row per (`TNR`, antibiotic) with a lab
  interpretation (`INT_ERG`/`ERG`: S/I/R/K=no criteria/U=unknown) and,
  for broth-microdilution results (Thermo Fisher Sensititre SLOMYCOI /
  SLOMYCO2 panels), an `MHK` (MIC) value. Rows where `ANTIBIOTIKA` encodes
  a fixed test concentration (e.g. `"Amikacin 4 mg/l"`) are single-point
  breakpoint tests and never carry an `MHK` value -- **this analysis only
  uses rows with a non-empty `MHK`**, per the intended scope.
- `data/imm/screening_map.csv` -- links `TNR` to sample/sequencing IDs and
  marks reference/control-strain rows (`LABEL`), including
  `M. kansasii (ATCC 12478)`. That control row has no `TNR`/MIC data of
  its own in this export, so it cannot serve as an in-dataset QC
  comparison -- see "Reference values" below. Its `MHK` column is a
  second, unrelated `TNR`-like id: for isolates that had their broth
  microdilution (MHK) test run under a different lab number than their
  primary `TNR`, that number is recorded there instead of in the `TNR`
  column.

**Joining the two files** (`src/kansasii_mic/loading.py`): a `mic.csv` row
is kept only if its `TNR` matches screening_map's `TNR` column or, failing
that, its `MHK` column -- screening_map defines which isolates are in
scope, so rows matching neither are dropped. On an `MHK`-column match, the
row's `TNR` is replaced with screening_map's canonical `TNR` for that
isolate. The resulting long-format table carries `NR`, `PROBENNUMMER`,
`TNR` and `MHK` from screening_map alongside the parsed MIC fields;
figures are labeled by `PROBENNUMMER` (a string), while `NR`/`TNR`/`MHK`
are nullable integers.

## Data-quality handling

- **Antibiotic name normalization**: `Sulfamethox.-Trimethop.` and
  `Trimethoprim/Sulfamethoxazol` are merged into one canonical
  `Sulfamethoxazole/Trimethoprim`; `Amikacin`, `Amikacin i.v.` and
  `Amikacin inhalativ` are merged into `Amikacin` (original label kept in
  `raw_antibiotic`).
- **`MHK` string parsing** (`src/kansasii_mic/parsing.py`) handles plain
  numbers, censored values (`<X`, `<=X`, `>X`), doubling-dilution ranges
  (`X-Y`), and fixed-ratio Sulfamethoxazole/Trimethoprim combo values
  (`X/Y`, including combo ranges), plus one confirmed data-entry typo
  (a combo range missing a slash), repaired generically rather than
  hardcoded.
- **Flat-profile artifact detection** (`src/kansasii_mic/qc.py`): at least
  one `TNR` (`2021311033`) has runs of many different antibiotics all
  reporting the exact same MIC value, which is not biologically plausible
  and looks like placeholder data. Any run of >= 8 distinct antibiotics
  within one `TNR`, in original file order, sharing an identical MIC point
  estimate is flagged and excluded from the main analysis (written to
  `excluded_qc.csv` with a reason) rather than hardcoding that one `TNR`.

## Point-estimate conventions

Doubling-dilution MICs are analyzed on the log2 scale.

- censored `<X` / `<=X` / `>X` -> boundary value `X` (documented
  simplification; the true value is only known to be beyond that bound)
- range `X-Y` -> geometric mean of the bounds (= arithmetic mean on log2)
- combo `X/Y` (Sulfamethoxazole/Trimethoprim, fixed 19:1 ratio) -> the
  first (trimethoprim numerator) component, matching how CLSI expresses
  that breakpoint (e.g. `2/38`)

## Reference values for *M. kansasii*

*M. kansasii* ATCC 12478 is the CLSI-designated quality-control strain for
the Sensititre SLOMYCO panels used to generate this data (primarily for
rifampin MIC QC), but it has no linked MIC results in this dataset, so
species-specific CLSI breakpoints are used instead
(`src/kansasii_mic/breakpoints.py`), summarized from CLSI M24/M62 and
van Ingen, "Antimycobacterial Susceptibility Testing of Nontuberculous
Mycobacteria" (PMC6760954):

| Antibiotic | S (<=) | I | R (>=) |
|---|---|---|---|
| Rifampicin | 1 | -- | 2 |
| Clarithromycin | 8 | 16 | 32 |
| Amikacin | 16 | 32 | 64 |
| Ciprofloxacin | 1 | 2 | 4 |
| Doxycyclin | 1 | 2-4 | 8 |
| Linezolid | 8 | 16 | 32 |
| Minocyclin | 1 | 2-4 | 8 |
| Moxifloxacin | 1 | 2 | 4 |
| Rifabutin | 2 | -- | 4 |
| Sulfamethoxazole/Trimethoprim | 2/38 | -- | 4/76 |

No CLSI *M. kansasii*-specific breakpoint exists for ethambutol,
isoniazid, streptomycin, ethionamide, bedaquiline or clofazimine --
consistent with this dataset's own `INT_ERG`/`ERG` marking those drugs
`K` ("keine Bewertung", no assessment criteria) for every isolate.

## Outlier detection

- **Per antibiotic** (`outliers_per_antibiotic.csv`): Tukey fences
  (median +/- 1.5 x IQR) on log2 MIC across all `TNR`s tested for that
  drug, flagging isolates as `more_resistant` or `more_susceptible` than
  the rest of the cohort. Robust to the small, skewed, discrete-dilution
  data here. Each row also carries the CLSI category (if defined) for
  cross-reference.
- **Per isolate, aggregated** (`tnr_resistance_ranking.csv`): mean robust
  z-score `(log2 MIC - median) / MAD` across all antibiotics that `TNR`
  was tested for, plus counts of CLSI-resistant results -- used to rank
  overall "most resistant" / "most susceptible" isolates.

Most antibiotics have MIC data for ~50-52 of the 54 usable isolates, but a
few are tested rarely (e.g. Clofazimin, n=2 in this dataset) -- Tukey
fences on that few points are uninformative, so treat outlier flags (or
their absence) for low-n antibiotics with caution; check
`n_antibiotics_tested`-style counts in the per-antibiotic figures/titles
before drawing conclusions.


## Running

```bash
module load miniforge3
conda env create -f environment.yml   # first time only
conda activate kansasii_mic

pytest tests/

python scripts/run_analysis.py \
  --mic-csv /shares/sander.imm.uzh/MM/kansasii/data/imm/mic.csv \
  --screening-map /shares/sander.imm.uzh/MM/kansasii/data/imm/screening_map.csv \
  --out-dir /shares/sander.imm.uzh/MM/kansasii/output/mic
```

## Outputs (`output/mic/`)

- `mic_parsed.csv` -- cleaned long-format MIC data (one row per isolate x
  antibiotic, with `NR`/`PROBENNUMMER`/`TNR`/`MHK` from screening_map)
  used for the analysis.
- `excluded_qc.csv` -- rows dropped as flat-profile data-quality artifacts.
- `parse_failures.csv` -- written only if any `MHK` value could not be
  parsed at all.
- `outliers_per_antibiotic.csv`, `tnr_resistance_ranking.csv` -- see above.
- `figures/<antibiotic>_distribution.png`, `figures/heatmap_tnr_antibiotic.png`,
  `figures/resistance_ranking.png`.
