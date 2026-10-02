# ==========================================================================
# example_single.metal  --  fixed-effect, inverse-variance meta-analysis
#                           for ONE phenotype/stratum.
# ==========================================================================

SCHEME   STDERR          # fixed-effect, inverse-variance weighted (beta + SE)
AVERAGEFREQ ON           # carry a weighted Effect_AF into the output
MINMAXFREQ  ON           # also report min/max freq across cohorts

# --- Define custom variables to track cases and controls ---
CUSTOMVARIABLE Ncases
CUSTOMVARIABLE Ncontrols

# --- column names in the QC'd files (from qc_munge.R) ---
MARKER    SNP_ID
ALLELE    Effect_allele Non-effect_allele
FREQLABEL Effect_AF
EFFECT    Beta
STDERR    SE
PVALUE    P-value

# --- Set the column labels ONCE globally ---
LABEL Ncases AS N_CAS
LABEL Ncontrols AS N_CON

# --- the cohorts to combine (add one PROCESS line per cohort) ---
PROCESS   BIOBANK1.tsv
PROCESS   BIOBANK2.tsv

OUTFILE   ATOPIC_DERM_EUR_ .tbl
ANALYZE HETEROGENEITY                 
QUIT
