# HTDV 0.2.0

Validation release. The framework now ships with the executable evidence
of its calibration claims.

## New: validation evidence shipped with the package

- `htdv_sim_summary`: aggregated output of the pre-registered factorial
  Monte Carlo study (1024 cells × 500 replications × 3 inferential
  layers, 31 hours wall-clock on 16 cores, sign-corrected coverage).
- `htdv_empirical_benchmarks`: three-dataset external validation against
  Stock-Watson (2007), Campbell-Shiller (1998), and the iid Welch
  baseline on US-Canada 10y yields.
- New vignette `HTDV-validation` consolidates both validations with
  the per-layer reading and the persistence-ladder cross-cut.

## New: simulation infrastructure

- `htdv_simstudy()`: end-to-end factorial Monte Carlo runner with
  cell-level parallelism via `parallel::mclapply` and atomic per-cell
  caching that survives worker failures and supports resume.
- `htdv_simstudy_summary()`: aggregate per-cell results.
- `htdv_simstudy_warnings()`: flag cells where the Bayesian
  diagnostic-pass rate falls below a user-set threshold (typically
  high-persistence × small-sample corners where the AR(1) likelihood
  approaches its limit of identification).
- New Stan model `simstudy_two_sample.stan` for two-sample AR(1)
  inference with shared autocorrelation.

## Bug fixes

- Sign convention in the simulation's HAR-Wald and stationary-bootstrap
  layers aligned with the Stan parameterization
  (`delta = E[X_2] - E[X_1]`). Earlier output had spuriously low HAR
  and bootstrap coverage for `delta > 0` cells; the
  `recover_coverage.R` helper of the companion paper repository
  regenerates the corrected aggregates from the raw RDS.
- Documentation references for Betancourt (2016) corrected to use
  `<doi:10.48550/arXiv.1604.00695>` (CRAN-mandated DOI form).
- Package metadata: GitHub URL points to `IsadoreNabi/HTDV`; ORCID,
  email, and BugReports URL aligned in DESCRIPTION and CITATION.

# HTDV 0.1.0

Initial public release.

- Hierarchical Bayesian fit via HMC (`htdv_fit`) with TAC, WSC, MPC, Whittle,
  and Composite likelihood backends.
- Berger-robust posterior envelope (`htdv_envelope`).
- Long-run variance with Andrews bandwidth (`htdv_lrv`).
- Fixed-b HAR Wald test (`htdv_fixedb`).
- Block bootstrap with Patton-Politis-White automatic block length
  (`htdv_boot`).
- Adaptive conformal inference for dependent data (`htdv_conformal`).
- ROPE decision, bridge-sampling Bayes factor, WAIC, leave-future-out CV,
  and predictive stacking (`htdv_rope`, `htdv_bf`, `htdv_waic_lfo`,
  `htdv_stack`).
- MCMC diagnostics and posterior-predictive checks on dependence statistics
  (`htdv_diagnostics`, `htdv_ppc`).
- Explicit numerical constants for TAC/WSC/MPC metric equivalence
  (`htdv_equivalence_constants`).
