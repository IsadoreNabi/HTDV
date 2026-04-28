# HTDV

Hypothesis Testing for Dependent Variables with Unbalanced Data.

`HTDV` provides a unified R toolkit for inference on dependent,
unbalanced data under strong-mixing conditions, combining hierarchical
Bayesian estimation via Hamiltonian Monte Carlo with frequentist and
distribution-free robustness anchors (fixed-b HAR, block bootstrap,
adaptive conformal).

## Validation card

The framework is shipped with two pre-registered validation studies,
both reproducible end-to-end and with their summary tables exposed as
package datasets. See `vignette("HTDV-validation")`.

- **Factorial Monte Carlo** (`htdv_sim_summary`). 1024-cell design
  crossing sample size, AR(1) coefficient, innovation tail, imbalance
  ratio and location shift; 500 replications per cell × 3 inferential
  layers; 31 hours of wall-clock on 16 cores. The Bayesian envelope
  holds nominal size (mean 0.0556, sd 0.013) and nominal coverage (mean
  0.944) across the entire grid; HAR and bootstrap inflate to
  empirical size 0.60 and coverage 0.29 in the worst corners under
  strong persistence. The asymptotic gap that motivates the framework
  is *visible in the data*.
- **External benchmarks** (`htdv_empirical_benchmarks`). Three public
  datasets compared against published references:
  - FRED-MD post-1984 CPI inflation against Stock and Watson (2007).
  - Shiller log-CAPE against Campbell and Shiller (1998).
  - US-Canada 10-year yield differential against the iid Welch baseline.

  All three layers reproduce all three references with `agreement` in
  every case. The 95% interval widths scale monotonically with the
  series persistence: at $\widehat\phi\approx 0.45$ Bayes is 0.81× HAR;
  at $\widehat\phi\approx 0.97$ it is 2.80× HAR; at near-unit-root
  ($\widehat\phi\approx 0.99$) it is 15.0× HAR. The framework's value
  is the *visibility* of this gradient.

```{r, eval = FALSE}
library(HTDV)
data(htdv_sim_summary)         # simulation summary, 3069 rows
data(htdv_empirical_benchmarks) # three-dataset external validation
vignette("HTDV-validation")    # full narrative
```

## Installation

```r
remotes::install_github("IsadoreNabi/HTDV")
```

`rstan` is required. Optional backends: `bridgesampling` (Bayes factors),
`loo` (WAIC / PSIS-LOO), `posterior` (draws utilities), `bayesplot`
(visualization), `readxl` (vignette).

## Core API

| Function | Purpose |
|----------|---------|
| `htdv_fit()` | Hierarchical Bayesian HMC fit. |
| `htdv_envelope()` | Berger-robust envelope across models. |
| `htdv_lrv()` | HAC long-run variance (Andrews bandwidth). |
| `htdv_fixedb()` | Fixed-bandwidth HAR Wald test. |
| `htdv_boot()` | Block bootstrap with automatic block length. |
| `htdv_conformal()` | Adaptive conformal inference. |
| `htdv_rope()` | ROPE-based posterior decision. |
| `htdv_bf()` | Bridge-sampling Bayes factor. |
| `htdv_waic_lfo()` | WAIC and leave-future-out CV. |
| `htdv_stack()` | Predictive stacking. |
| `htdv_diagnostics()` | MCMC diagnostics. |
| `htdv_ppc()` | Posterior-predictive checks on dependence statistics. |
| `htdv_equivalence_constants()` | Explicit TAC/WSC/MPC constants. |
| `htdv_simstudy()` | Factorial Monte Carlo study (Section 12-bis). |
| `htdv_simstudy_summary()` | Aggregate per-cell results. |
| `htdv_simstudy_warnings()` | Flag cells in the limit-of-identification zone. |

See `vignette("HTDV-intro")` for a walkthrough,
`vignette("HTDV-validation")` for the full validation report.

## Citation

Please cite both the package and the companion paper. Run
`citation("HTDV")` for the current BibTeX entries.

## License

MIT.
