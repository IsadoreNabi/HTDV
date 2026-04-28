# HTDV 0.1.0 submission

## Test environments

- local: Fedora 43, R 4.4.x
- win-builder (release, devel)
- macOS release via R-hub

## R CMD check results

0 errors | 0 warnings | 0 notes (expected on submission).

## Stan compilation

The package defers Stan model compilation to first invocation rather than
installing pre-compiled binaries. All Stan programs ship in `inst/stan/` and
are compiled via `rstan::stan_model()` on demand; compiled objects are
cached for the session in a package-private environment.

## Suggested-package usage

All Suggests (`bridgesampling`, `loo`, `posterior`, `bayesplot`) are accessed
through `requireNamespace(..., quietly = TRUE)` guards; the package functions
that require them emit a clear message if the suggested package is not
installed and never attempt to install packages on behalf of the user.

## Reverse dependencies

None on initial release.
