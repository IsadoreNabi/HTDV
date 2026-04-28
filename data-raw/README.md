# `data-raw/` — reproducibility scripts

This directory contains the executable scripts that produced the
package datasets (`htdv_sim_summary`, `htdv_empirical_benchmarks`)
and supplied the empirical sections of the companion paper. They are
shipped in the GitHub source tree but excluded from the CRAN tarball
via `.Rbuildignore`.

| Script | Purpose |
|--------|---------|
| `run_simstudy.R` | Full-factorial 1024-cell Monte Carlo (Section 12-bis of the paper). Wall-clock ~31h on 16 cores. |
| `recover_coverage.R` | Post-hoc sign correction of HAR/bootstrap coverage in the simulation output. ~1 min. |
| `run_benchmarks.R` | Three-dataset external validation (Section 12-ter): FRED-MD inflation, Shiller log-CAPE, US-Canada 10y yield differential. ~10-20 min. |
| `build_package_data.R` | Converts the two summary CSVs into compressed `.rda` artifacts in `data/`. |

## Required data files

The empirical benchmarks read three public files. Download manually and
save them in this directory (paths are relative to `data-raw/`):

| File | Source URL |
|------|------------|
| `2026-03-MD.csv` | https://www.stlouisfed.org/research/economists/mccracken/fred-databases |
| `ie_data.xls` | https://shillerdata.com/ |
| `us_canada_10y.csv` | https://fred.stlouisfed.org/graph/fredgraph.csv?id=GS10,IRLTLT01CAM156N |

The data files are not committed to the repository (see `.gitignore`).

## Reproducing from scratch

```bash
# 1. Run the simulation (long: ~31h on 16 cores)
Rscript run_simstudy.R

# 2. Apply the post-hoc sign correction
Rscript recover_coverage.R

# 3. Run the three empirical benchmarks (~10-20 min)
Rscript run_benchmarks.R

# 4. Convert summary CSVs to package datasets
Rscript build_package_data.R

# 5. Reinstall the package with the refreshed datasets
R -e "devtools::install('.', build_vignettes = TRUE, upgrade = FALSE)"
```
