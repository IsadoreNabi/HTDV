# =====================================================================
# build_package_data.R
# -----------------------------------------------------------------
# Converts the two summary CSVs produced by the simulation and the
# empirical-benchmark runs into compressed .rda artifacts inside
# HTDV/data/ so that they ship with the package and become available
# via data() after install.
#
# Run from R with:
#   source("build_package_data.R")
# or from the shell with:
#   Rscript build_package_data.R
# =====================================================================

.resolve_script_dir <- function() {
  frame_ofile <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(frame_ofile) && nzchar(frame_ofile))
    return(dirname(normalizePath(frame_ofile)))
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m) == 1L)
    return(dirname(normalizePath(sub("^--file=", "", m))))
  getwd()
}
out_dir <- .resolve_script_dir()

sim_csv <- file.path(out_dir, "simstudy_summary_v2.csv")
emp_csv <- file.path(out_dir, "benchmarks_E1_E2_E3.csv")
data_dir <- file.path(out_dir, "HTDV", "data")

if (!file.exists(sim_csv)) stop("Missing: ", sim_csv)
if (!file.exists(emp_csv)) stop("Missing: ", emp_csv)
if (!dir.exists(data_dir))
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading: ", sim_csv)
htdv_sim_summary <- read.csv(sim_csv, header = TRUE,
                             stringsAsFactors = FALSE)
message(sprintf("  %d rows x %d columns", nrow(htdv_sim_summary),
                ncol(htdv_sim_summary)))

message("Reading: ", emp_csv)
htdv_empirical_benchmarks <- read.csv(emp_csv, header = TRUE,
                                      stringsAsFactors = FALSE)
message(sprintf("  %d rows x %d columns", nrow(htdv_empirical_benchmarks),
                ncol(htdv_empirical_benchmarks)))

save(htdv_sim_summary,
     file = file.path(data_dir, "htdv_sim_summary.rda"),
     compress = "xz")
save(htdv_empirical_benchmarks,
     file = file.path(data_dir, "htdv_empirical_benchmarks.rda"),
     compress = "xz")

message("\nWrote:")
message("  ", file.path(data_dir, "htdv_sim_summary.rda"),
        sprintf(" (%.1f KB)",
                file.info(file.path(data_dir, "htdv_sim_summary.rda"))$size / 1024))
message("  ", file.path(data_dir, "htdv_empirical_benchmarks.rda"),
        sprintf(" (%.1f KB)",
                file.info(file.path(data_dir, "htdv_empirical_benchmarks.rda"))$size / 1024))

message("\nReinstall the package and the new datasets and vignette ship with it:")
message("  devtools::document('HTDV')")
message("  devtools::install('HTDV', build_vignettes = TRUE, upgrade = 'never')")
