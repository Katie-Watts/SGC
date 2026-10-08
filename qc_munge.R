#!/usr/bin/env Rscript
###############################################################################
# qc_munge.R  --  QC / harmonise cohort-level GWAS summary statistics with
#                 MungeSumstats, then emit the consortium-standard columns.
#
# The process (per cohort file)
#   1. Runs MungeSumstats::format_sumstats() with the SGC QC filters
#      (INFO >= 0.3, FRQ >= 0.005; MAF >= 0.005 is applied in step 2). Column headers are auto-detected from
#      MungeSumstats' default map (effect_AF is mapped to FRQ explicitly as not a default option), so
#      cohort files don't need pre-renaming. Each cohort's genome build is
#      inferred and lifted to GRCh38; missing RSIDs are filled from dbSNP 155
#      where found. Variants not in dbSNP are kept; variants that fail liftover are dropped.
#   2. Reads the munged file back and rewrites it with standard column names
#      plus a leading SNP_ID column
#      (Chromosome:Position:Allele1:Allele2, alleles sorted alphabetically) that METAL uses
#      as the marker key in the next step. N_CAS / N_CON are passed through
#      when the cohort supplies them, otherwise written as NA.
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

# QC filters to impose(SGC standard)
INFO_FILTER <- 0.3
FRQ_FILTER  <- 0.005

# Genome build handling: infer each cohort's build, then lift everything to
# GRCh38 so all cohorts share coordinates before meta-analysis (as a backup check as Drew's code should handle this)
REF_GENOME_DEFAULT <- NULL       # NULL = let MungeSumstats infer per file
CONVERT_REF_TO     <- "GRCh38"   # target build for all outputs
DBSNP_BUILD        <- 155        # dbSNP reference for SNP mapping

# Keep variants not found in dbSNP / the reference; RSIDs are only added
# where found (and not already present).
DROP_NON_DBSNP     <- FALSE

# ============================================================================

STRATA <- c("ALL","EUR","AFR","AMR","EAS","SAS","MALE","FEMALE")

# CLI args
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

# input files
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

# header map
# NB: in MungeSumstats A2 is the EFFECT allele (Can rename to anything we want though.)
mapping_file <- MungeSumstats::sumstatsColHeaders

# remapping as EFFECT_AF not a default mapping (effect-allele freq)
mapping_file <- unique(rbind(mapping_file,
                             data.frame(Uncorrected = "EFFECT_AF", Corrected = "FRQ")))

# remapping for SE
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
# Need to check NA codings and capitalisation
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
  c_cas <- pick(dt, "N_CAS")   # optional: written as NA if the cohort didn't supply it
  c_con <- pick(dt, "N_CON")
  c_snp <- pick(dt, "SNP")     # rsID from dbSNP where found (MungeSumstats SNP column)

  # MAF filter: FRQ is the effect-allele freq, so fold it to catch rare variants
  # whichever allele is the effect allele. Rows with no FRQ are kept.
  if (!is.na(c_frq)) {
    frq  <- as.numeric(dt[[c_frq]])
    keep <- is.na(frq) | pmin(frq, 1 - frq) >= FRQ_FILTER
    if (any(!keep))
      message(sprintf("    [MAF] dropped %s variant(s) with MAF < %s",
                      format(sum(!keep), big.mark = ","), FRQ_FILTER))
    dt <- dt[keep]
  }

  chr <- sub("^chr", "", as.character(dt$CHR), ignore.case = TRUE)
  ea  <- toupper(dt$A2)   # A2 = effect allele (MungeSumstats convention)
  oa  <- toupper(dt$A1)   # A1 = non-effect / reference allele

  # Orientation-independent marker key: alleles sorted alphabetically so the
  # same variant gets the same SNP_ID whichever allele a cohort calls "effect".
  # METAL aligns betas itself using Effect_allele / Non-effect_allele.
  a_lo <- pmin(ea, oa)
  a_hi <- pmax(ea, oa)

  #Rename ouptut columns as needed / whatever we want but will then need updating in METAL script if changed
  out <- data.table::data.table( 
    SNP_ID              = paste(chr, dt$BP, a_lo, a_hi, sep = ":"),
    RSID                = if (!is.na(c_snp)) as.character(dt[[c_snp]]) else NA_character_,
    Chromosome          = chr,
    Position            = dt$BP,
    Effect_allele       = ea,
    `Non-effect_allele` = oa,
    Beta                = dt$BETA,     # relative to A2
    SE                  = dt$SE,
    `P-value`           = dt$P,
    Effect_AF           = if (!is.na(c_frq)) as.numeric(dt[[c_frq]]) else NA_real_,  # A2 freq
    Imp_Quality         = if (!is.na(c_inf)) dt[[c_inf]] else NA_real_,
    N_CAS               = if (!is.na(c_cas)) as.numeric(dt[[c_cas]]) else NA_real_,
    N_CON               = if (!is.na(c_con)) as.numeric(dt[[c_con]]) else NA_real_
  )
  # Write to a temp file and rename only once complete, so a crash mid-write
  # never leaves a truncated output that the resume check would then skip.
  tmp <- paste0(out_path, ".tmp")
  data.table::fwrite(out, tmp, sep = "\t", quote = FALSE, na = "NA",
                     compress = "gzip")
  if (!file.rename(tmp, out_path))
    stop("could not rename ", tmp, " to ", out_path)
  nrow(out)
}

# ------------------------------------------------------------------ main loop
n_ok <- 0L; n_skip <- 0L; n_fail <- 0L; checked <- FALSE
for (f in files) {
  base <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", basename(f))
  out_path <- file.path(OUT, paste0(base, ".tsv.gz"))
  if (file.exists(out_path) && !FORCE) {
    message(sprintf("[skip] %s (output exists)", base)); n_skip <- n_skip + 1L; next
  }

  message(sprintf("[munge] %s", base))
  munged_path <- file.path(munged_dir, paste0(base, ".munged.tsv.gz"))
  res <- tryCatch(
    MungeSumstats::format_sumstats(
      path                 = f,
      save_path            = munged_path,
      ref_genome           = REFG,
      convert_ref_genome   = CONVERT_REF_TO,
      dbSNP                = DBSNP_BUILD,
      on_ref_genome        = DROP_NON_DBSNP,   # FALSE = keep SNPs not found in dbSNP
      INFO_filter          = INFO_FILTER,
      FRQ_filter           = FRQ_FILTER,
      bi_allelic_filter    = TRUE,
      allele_flip_check    = TRUE,
      N_dropNA             = FALSE,
      snp_ids_are_rs_ids   = FALSE,
      mapping_file         = mapping_file,
      nThread              = THREADS,
      log_folder           = file.path(log_root, base),
      log_folder_ind       = TRUE,
      force_new            = FORCE,
      return_data          = FALSE
    ),
    error = function(e) { message("    [FAIL munge] ", conditionMessage(e)); NULL }
  )
  if (is.null(res)) { n_fail <- n_fail + 1L; next }

  # format_sumstats returns the save path (character) or, with log_folder_ind,
  # a list whose $sumstats is that path.
  mp <- if (is.list(res)) res$sumstats else res
  if (is.null(mp) || !file.exists(mp)) {
    message("    [FAIL] munged file not found"); n_fail <- n_fail + 1L; next
  }

  n <- tryCatch(standardise(mp, out_path),
                error = function(e) { message("    [FAIL standardise] ",
                                              conditionMessage(e)); NA })
  if (is.na(n)) { n_fail <- n_fail + 1L; next }
  message(sprintf("    [ok] %s  (%s variants)", basename(out_path),
                  format(n, big.mark = ",")))
  n_ok <- n_ok + 1L
}

message(sprintf("\nDONE.  %d munged, %d skipped, %d failed.  Output -> %s/",
                n_ok, n_skip, n_fail, OUT))
if (n_fail > 0) message("See per-file logs under ", log_root, "/ for failures.")
