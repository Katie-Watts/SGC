# ==========================================================================
# example_single.metal  --  fixed-effect, inverse-variance meta-analysis
#                           for ONE phenotype/stratum (edit for your run).
#
# Run:
#   # METAL cannot read gzip -- decompress the QC'd inputs first:
#   gzip -dc qc/ATOPIC_DERM_EUR_biobankX.tsv.gz > biobankX.tsv
#   gzip -dc qc/ATOPIC_DERM_EUR_cohortY.tsv.gz   > cohortY.tsv
#   metal < example_single.metal
# ==========================================================================

SCHEME   STDERR          # fixed-effect, inverse-variance weighted (beta + SE)
AVERAGEFREQ ON           # carry a weighted Effect_AF into the output
MINMAXFREQ  ON           # also report min/max freq across cohorts

# --- column names in the QC'd files (from qc_munge.R) ---
MARKER    SNP_ID
ALLELE    Effect_allele Non-effect_allele
FREQLABEL Effect_AF
EFFECT    Beta
STDERR    SE
PVALUE    P-value

# --- the cohorts to combine (add one PROCESS line per cohort) ---
PROCESS   biobankX.tsv
PROCESS   cohortY.tsv

OUTFILE   ATOPIC_DERM_EUR_ .tbl
ANALYZE                  # use  ANALYZE HETEROGENEITY  for I^2 / Q stats
QUIT
