#!/usr/bin/env Rscript
###############################################################################
# qc_munge.R  --  QC / harmonise cohort-level GWAS summary statistics with
#                 MungeSumstats, then emit the consortium-standard columns.
#
# The process (per cohort file)
#   1. Runs MungeSumstats::format_sumstats() with the SGC QC filters
#      (INFO >= 0.3, FRQ >= 0.005). Column headers are auto-detected from
#      MungeSumstats' default map (effect_AF is mapped to FRQ explicitly as not a default option), so
#      cohort files don't need pre-renaming. Each cohort's genome build is
#      inferred and lifted to GRCh38; missing RSIDs are filled from dbSNP 155
#      (variants not in dbSNP are still kept).
#   2. Reads the munged file back and rewrites it with standard column names
#      plus a leading SNP_ID column
#      (Chromosome:Position:Effect_allele:Non-effect_allele) that METAL uses
#      as the marker key in the next step.
#   3. Writes <out_dir>/<same basename>.tsv.gz (tab-separated, gzipped).
#      Drew - Rename files if needed.
#
# RESUMABLE: a cohort whose output already exists is skipped (unless --force).
#
# ONE-OFF SETUP (installs MungeSumstats + references) -- run once, separately:
#   Rscript qc_setup.R
#
# USAGE
#   Rscript qc_munge.R --in cohort_files --out qc [options]
#   Rscript qc_munge.R --in 'cohort_files/*.txt.gz' --out qc --threads 4
#
# OPTIONS
#   --in PATH           directory, single file, or quoted glob (default: cohort_files)
#   --out DIR           output directory (default: qc)
#   --threads N         threads for MungeSumstats (default: 1)
#   --mapping FILE      extra header map (.csv/.tsv/.xlsx; Uncorrected,Corrected)
#   --ref-genome BUILD  force input build (GRCh37/GRCh38); default: infer per file
#   --force             re-run cohorts whose output already exists
#   --check N           (default 0 = off)
#   -h, --help          show this header
###############################################################################

# QC filters (SGC standard)
INFO_FILTER <- 0.3
FRQ_FILTER  <- 0.005

# Genome build handling: infer each cohort's build, then lift everything to
# GRCh38 so all cohorts share coordinates before meta-analysis (backup only;
# cohorts were asked to supply GRCh38).
REF_GENOME_DEFAULT <- NULL       # NULL = let MungeSumstats infer per file
CONVERT_REF_TO     <- "GRCh38"   # target build for all outputs
DBSNP_BUILD        <- 155        # dbSNP reference for SNP mapping

# Keep variants not found in dbSNP / the reference; RSIDs are only added
# where found (and not already present).
DROP_NON_DBSNP     <- FALSE

# ============================================================================

STRATA <- c("ALL","EUR","AFR","AMR","EAS","SAS","MALE","FEMALE")

# ------------------------------------------------------------------ CLI args
args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(sprintf("%s needs a value", flag))
  args[i + 1]
}
has_flag <- function(flag) flag %in% args

if (has_flag("--help") || has_flag("-h")) {
  lines <- readLines(sub("--file=", "",
      grep("--file=", commandArgs(FALSE), value = TRUE)))
  hdr_end <- which(!grepl("^#", lines))[1] - 1   # print the leading comment block only
  cat(lines[2:hdr_end], sep = "\n")
  quit(status = 0)
}

IN      <- get_opt("--in",  default = "cohort_files")
OUT     <- get_opt("--out", default = "qc")
THREADS <- as.integer(get_opt("--threads", default = "1"))
MAPPING <- get_opt("--mapping", default = NULL)   # optional extra header map
REFG    <- get_opt("--ref-genome", default = REF_GENOME_DEFAULT)
FORCE   <- has_flag("--force")
CHECK_N <- as.integer(get_opt("--check", default = "0"))

# required packages (installed once by qc_setup.R)
for (p in c("MungeSumstats", "data.table")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("Package '%s' is not installed. Run the setup first:  Rscript qc_setup.R", p))
}
suppressWarnings(suppressMessages(library(data.table)))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
munged_dir <- file.path(OUT, "_munged"); dir.create(munged_dir, showWarnings = FALSE)
log_root   <- file.path(OUT, "_logs");   dir.create(log_root,   showWarnings = FALSE)

# ------------------------------------------------------------------ input files
if (grepl("[*?]", IN)) {
  files <- Sys.glob(IN)
} else if (dir.exists(IN)) {
  files <- list.files(IN, pattern = "\\.(tsv|txt|csv)(\\.gz)?$", full.names = TRUE)
} else if (file.exists(IN)) {
  files <- IN
} else {
  stop("No input found at: ", IN)
}
files <- files[!grepl("(^|/)\\._", files)]        # skip AppleDouble junk
if (!length(files)) stop("No cohort files matched: ", IN)
message(sprintf("Found %d cohort file(s).", length(files)))

# ------------------------------------------------------------------ header map
# MungeSumstats' default map handles the SGC headers (effect_allele -> A2,
# non_effect_allele -> A1). NB: in MungeSumstats A2 is the EFFECT allele, and
# BETA / FRQ are relative to A2 after alignment to the reference.
mapping_file <- MungeSumstats::sumstatsColHeaders
# make sure the SGC frequency header is recognised as FRQ (effect-allele freq)
mapping_file <- unique(rbind(mapping_file,
                             data.frame(Uncorrected = "EFFECT_AF", Corrected = "FRQ")))

# optional extra mappings (--mapping) for oddities the defaults miss,
# e.g. a "standard error" column with a space. Two columns: Uncorrected, Corrected.
if (!is.null(MAPPING)) {
  extra <- if (tolower(tools::file_ext(MAPPING)) %in% c("xlsx","xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("readxl needed to read ", MAPPING)
    as.data.frame(readxl::read_excel(MAPPING))
  } else {
    as.data.frame(data.table::fread(MAPPING))
  }
  names(extra)[1:2] <- c("Uncorrected", "Corrected")
  extra$Uncorrected <- toupper(extra$Uncorrected)
  mapping_file <- unique(rbind(mapping_file, extra[, c("Uncorrected","Corrected")]))
  message(sprintf("Loaded %d extra header mapping(s) from %s", nrow(extra), MAPPING))
}

# ------------------------------------------------------------------ standardise
# Munged files always use MungeSumstats' standard names, so only a few
# fallbacks are needed.
pick <- function(dt, candidates) {
  hit <- intersect(candidates, names(dt))
  if (length(hit)) hit[1] else NA_character_
}

standardise <- function(munged_path, out_path) {
  dt <- data.table::fread(munged_path, sep = "\t", header = TRUE,
                          na.strings = c("NA","na","NaN",""))

  need <- c("CHR","BP","A1","A2","BETA","SE","P")
  miss <- setdiff(need, names(dt))
  if (length(miss))
    stop("munged file missing required column(s): ", paste(miss, collapse = ", "))

  c_frq <- pick(dt, "FRQ")
  c_inf <- pick(dt, "INFO")

  chr <- sub("^chr", "", as.character(dt$CHR), ignore.case = TRUE)
  ea  <- toupper(dt$A2)   # A2 = effect allele (MungeSumstats convention)
  oa  <- toupper(dt$A1)   # A1 = non-effect / reference allele

  out <- data.table::data.table(
    SNP_ID              = paste(chr, dt$BP, ea, oa, sep = ":"),
    Chromosome          = chr,
    Position            = dt$BP,
    Effect_allele       = ea,
    `Non-effect_allele` = oa,
    Beta                = dt$BETA,     # relative to A2
    SE                  = dt$SE,
    `P-value`           = dt$P,
    Effect_AF           = if (!is.na(c_frq)) as.numeric(dt[[c_frq]]) else NA_real_,  # A2 freq
    Imp_Quality         = if (!is.na(c_inf)) dt[[c_inf]] else NA_real_
  )
  data.table::fwrite(out, out_path, sep = "\t", quote = FALSE, na = "NA",
                     compress = "gzip")
  nrow(out)
}
