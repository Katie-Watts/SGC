#!/usr/bin/env Rscript
###############################################################################
# qc_setup.R  --  ONE-OFF setup for qc_munge.R.
#
# Installs MungeSumstats, data.table, and the dbSNP / genome reference packages
# it needs to map SNPs and harmonise builds. These are large (several GB) and
# download once; after this, qc_munge.R runs offline against them.
#
# Run once per machine (or R library):
#   Rscript qc_setup.R
#
# Safe to re-run: packages already present are reported "ok" and not reinstalled
# (pass --force to reinstall them anyway).
###############################################################################

args  <- commandArgs(trailingOnly = TRUE)
FORCE <- "--force" %in% args

pkgs <- c(
  "MungeSumstats", "data.table",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38",
  "BSgenome.Hsapiens.NCBI.GRCh38",
  "SNPlocs.Hsapiens.dbSNP155.GRCh37",
  "BSgenome.Hsapiens.1000genomes.hs37d5"
)

message(">> Installing MungeSumstats and its references (one-off, large).")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

fail <- character(0)
for (p in pkgs) {
  if (!FORCE && requireNamespace(p, quietly = TRUE)) {
    message("   ok        ", p)
    next
  }
  message("   installing ", p)
  ok <- tryCatch({ BiocManager::install(p, update = FALSE, ask = FALSE); TRUE },
                 error = function(e) { message("      [FAIL] ", conditionMessage(e)); FALSE })
  if (!ok || !requireNamespace(p, quietly = TRUE)) fail <- c(fail, p)
}

if (length(fail)) {
  message("\n>> Setup INCOMPLETE. These packages did not install:")
  for (p in fail) message("     - ", p)
  message(">> Fix the errors above and re-run:  Rscript qc_setup.R")
  quit(status = 1)
}

message("\n>> Setup complete. You can now run:  Rscript qc_munge.R --in cohort_files --out qc")
