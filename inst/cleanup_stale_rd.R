# Run this ONCE before devtools::document() if upgrading from a prior version.
# Removes stale auto-generated .Rd files so the fresh correct ones are written.
# Usage (from RStudio with package root as working dir):
#   source("inst/cleanup_stale_rd.R")

man_dir <- "man"
if (!dir.exists(man_dir)) {
  message("man/ directory not found — nothing to clean.")
} else {
  old_rds <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  if (length(old_rds) == 0) {
    message("No stale .Rd files found.")
  } else {
    message("Removing ", length(old_rds), " stale .Rd file(s):")
    message(paste(" -", basename(old_rds), collapse = "\n"))
    invisible(file.remove(old_rds))
    message("\nDone. Now run devtools::document() then devtools::install().")
  }
}
