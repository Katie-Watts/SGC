# ==========================================================================
# example_single.metal  --  fixed-effect, inverse-variance meta-analysis
#                           for ONE phenotype/stratum.
#
# Run:
#   # METAL cannot read gzip -- decompress the QC'd inputs first:
#   gzip -dc qc/ATOPIC_DERM_EUR_BIOBANK1.tsv.gz > BIOBANK1.tsv
#   gzip -dc qc/ATOPIC_DERM_EUR_BIOBANK2.tsv.gz   > BIOBANK2.tsv
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
PROCESS   BIOBANK1.tsv
PROCESS   BIOBANK2.tsv

OUTFILE   ATOPIC_DERM_EUR_ .tbl
ANALYZE                  # HETEROGENEITY stats
QUIT
