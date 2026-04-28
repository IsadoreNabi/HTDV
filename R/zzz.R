.htdv_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  assign("stan_cache", list(), envir = .htdv_env)
  invisible(NULL)
}

.onUnload <- function(libpath) {
  if (exists("stan_cache", envir = .htdv_env, inherits = FALSE)) {
    rm("stan_cache", envir = .htdv_env)
  }
  invisible(NULL)
}
