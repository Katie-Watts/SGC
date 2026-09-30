#!/usr/bin/env Rscript
###############################################################################
# qc_munge.R  --  QC / harmonise cohort-level GWAS summary statistics with
#                 MungeSumstats, then emit the consortium-standard columns.
#
# The process (per cohort file)
#   1. Runs MungeSumstats::format_sumstats() with the SGC QC filters
#      (INFO >= 0.3, FRQ >= 0.005), harmonising to a single
#      allele convention. MungeSumstats auto-detects each cohort's (differing)
#      column headers, so we do not have to pre-rename them.
#   2. Reads the munged file back and rewrites it with exact standard
#      column names every SGC file must have:
#         Chromosome  Position  Effect_allele  Non-effect_allele
#         Beta  SE  P-value  Effect_AF  Imp_Quality
#      plus a leading SNP_ID column (Chromosome:Position:Effect_allele:
#      Non-effect_allele) that METAL uses as the marker key in the next step.
#   3. Writes <out_dir>/<same basename>.tsv.gz  (tab-separated, gzipped).
#      Keeping the input basename is what lets run_metal.py group cohorts of
#      the same PHENOTYPE_STRATUM together.
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
# ---------------------------------------------------------------------------
# ALLELE HANDLING
#   The cohort files already label the effect allele and its frequency
#   explicitly (columns  effect_allele  and  effect_AF ). We recode those to
#   MungeSumstats' own names on the way IN -- effect_allele -> A1, the other
#   allele -> A2, effect_AF -> FRQ -- so A1 is the effect allele and FRQ is the
#   effect-allele frequency throughout. MungeSumstats then keeps BETA and FRQ
#   oriented to A1 as it harmonises, and the output is taken AS-IS: no A1/A2
#   guessing, no sign flip, no 1-FRQ. Adjust the header names below if your
#   cohorts spell these columns differently.
###############################################################################

suppressWarnings(suppressMessages({
  ok_dt <- requireNamespace("data.table", quietly = TRUE)
}))

# ===================== INPUT HEADER NAMES (recode -> MSS) ===================
# How the cohort files spell the effect allele, the other allele, and the
# effect-allele frequency. These are recoded to MungeSumstats' A1 / A2 / FRQ
# so A1 = effect allele and FRQ = effect-allele frequency in the output.
# Case-insensitive; add extra spellings if cohorts vary.
EFFECT_ALLELE_NAMES <- c("effect_allele", "EA", "ALLELE1", "A1")
OTHER_ALLELE_NAMES  <- c("non_effect_allele", "other_allele", "NEA",
                         "ALLELE0", "ALLELE2", "A2")
EFFECT_FREQ_NAMES   <- c("effect_AF", "effect_af", "EAF", "eaf",
                         "effect_allele_frequency", "A1FREQ")
# ============================================================================

# QC filters (SGC standard)
INFO_FILTER <- 0.3
FRQ_FILTER  <- 0.005

# Genome build handling: infer each cohort's build, then lift everything to
# GRCh38 so all cohorts share coordinates before meta-analysis.
REF_GENOME_DEFAULT     <- NULL       # NULL = let MungeSumstats infer per file
CONVERT_REF_TO         <- "GRCh38"   # target build for all outputs
DBSNP_BUILD            <- 155        # dbSNP reference for SNP mapping
# ============================================================================

STRATA <- c("ALL","EUR","AFR","AMR","EAS","SAS","MALE","FEMALE")  # for --check only

# ------------------------------------------------------------------ arg parse
args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(flag, default = NULL, has_val = TRUE) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (!has_val) return(TRUE)
  if (i == length(args)) stop(sprintf("%s needs a value", flag))
  args[i + 1]
}
has_flag <- function(flag) flag %in% args

if (has_flag("--help") || has_flag("-h")) {
  cat(readLines(sub("--file=", "",
      grep("--file=", commandArgs(FALSE), value = TRUE)))[2:44], sep = "\n")
  quit(status = 0)
}

IN      <- get_opt("--in",  default = "cohort_files")
OUT     <- get_opt("--out", default = "qc")
THREADS <- as.integer(get_opt("--threads", default = "1"))
MAPPING <- get_opt("--mapping", default = NULL)   # optional custom header map
REFG    <- get_opt("--ref-genome", default = REF_GENOME_DEFAULT)
FORCE   <- has_flag("--force")
CHECK_N <- as.integer(get_opt("--check", default = "0"))

# require the heavy packages (installed once by qc_setup.R)
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
  files <- list.files(IN, pattern = "\\.(tsv|txt|csv)(\\.gz)?$",
                      full.names = TRUE)
} else if (file.exists(IN)) {
  files <- IN
} else {
  stop("No input found at: ", IN)
}
files <- files[!grepl("(^|/)\\._", files)]        # skip AppleDouble junk
if (!length(files)) stop("No cohort files matched: ", IN)
message(sprintf("Found %d cohort file(s).", length(files)))

# Header mapping. Start from MungeSumstats' default map, then RECODE the
# cohort's effect-allele / other-allele / effect-freq spellings to A1 / A2 /
# FRQ so A1 is the effect allele and FRQ is the effect-allele frequency.
mapping_file <- MungeSumstats::sumstatsColHeaders
mapping_file$Uncorrected <- as.character(mapping_file$Uncorrected)
mapping_file$Corrected   <- as.character(mapping_file$Corrected)
recode <- data.frame(
  Uncorrected = toupper(c(EFFECT_ALLELE_NAMES, OTHER_ALLELE_NAMES, EFFECT_FREQ_NAMES)),
  Corrected   = c(rep("A1",  length(EFFECT_ALLELE_NAMES)),
                  rep("A2",  length(OTHER_ALLELE_NAMES)),
                  rep("FRQ", length(EFFECT_FREQ_NAMES))),
  stringsAsFactors = FALSE)
# our recodes take precedence over the defaults for the same source header
mapping_file <- mapping_file[!(toupper(mapping_file$Uncorrected) %in% recode$Uncorrected), ]
mapping_file <- unique(rbind(recode, mapping_file[, c("Uncorrected","Corrected")]))

# optional EXTRA custom header mapping (--mapping), for oddities the auto-mapper
# misses, e.g. a "standard error" column WITH a space -- add a row  standard error -> SE .
if (!is.null(MAPPING)) {
  ext <- tools::file_ext(MAPPING)
  extra <- if (tolower(ext) %in% c("xlsx","xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("readxl needed to read ", MAPPING)
    as.data.frame(readxl::read_excel(MAPPING))
  } else {
    as.data.frame(data.table::fread(MAPPING))
  }
  names(extra)[1:2] <- c("Uncorrected", "Corrected")
  extra$Uncorrected <- toupper(extra$Uncorrected)
  mapping_file <- unique(rbind(mapping_file, extra[, c("Uncorrected","Corrected")]))
  message(sprintf("Loaded %d extra header mapping(s) from %s",
                  nrow(extra), MAPPING))
}

# ------------------------------------------------------------------ helpers
pick <- function(dt, candidates) {
  # first present column name (case-insensitive) from candidates, else NA
  for (c in candidates) {
    hit <- names(dt)[toupper(names(dt)) == toupper(c)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

standardise <- function(munged_path, out_path) {
  dt <- data.table::fread(munged_path, sep = "\t", header = TRUE,
                          na.strings = c("NA","na","NaN",""))
  up <- toupper(names(dt))

  c_chr <- pick(dt, c("CHR","CHROMOSOME"))
  c_bp  <- pick(dt, c("BP","POS","POSITION","BP_GRCH38"))
  c_a1  <- pick(dt, c("A1"))
  c_a2  <- pick(dt, c("A2"))
  c_b   <- pick(dt, c("BETA"))
  c_se  <- pick(dt, c("SE","STANDARD_ERROR"))
  c_p   <- pick(dt, c("P","PVAL","PVALUE","P_VALUE"))
  c_frq <- pick(dt, c("FRQ","FREQ","EAF","MAF","A2FREQ","A1FREQ"))
  c_inf <- pick(dt, c("INFO","IMPINFO","IMP_QUALITY","RSQ","R2"))

  need <- c(CHR=c_chr, BP=c_bp, A1=c_a1, A2=c_a2, BETA=c_b, SE=c_se, P=c_p)
  if (any(is.na(need)))
    stop("munged file missing required column(s): ",
         paste(names(need)[is.na(need)], collapse = ", "))

  # A1 = effect allele, A2 = other allele (recoded on input), so read directly.
  eff <- c_a1
  oth <- c_a2

  # Beta is always present and MungeSumstats keeps it oriented to A1.
  beta <- dt[[c_b]]

  # Effect_AF = FRQ (already the A1 / effect-allele frequency).
  eaf <- if (!is.na(c_frq)) as.numeric(dt[[c_frq]]) else NA_real_

  chr <- sub("^chr", "", as.character(dt[[c_chr]]), ignore.case = TRUE)
  ea  <- toupper(as.character(dt[[eff]]))
  oa  <- toupper(as.character(dt[[oth]]))

  out <- data.table::data.table(
    SNP_ID            = paste(chr, dt[[c_bp]], ea, oa, sep = ":"),
    Chromosome        = chr,
    Position          = dt[[c_bp]],
    Effect_allele     = ea,
    `Non-effect_allele` = oa,
    Beta              = beta,
    SE                = dt[[c_se]],
    `P-value`         = dt[[c_p]],
    Effect_AF         = eaf,
    Imp_Quality       = if (!is.na(c_inf)) dt[[c_inf]] else NA_real_
  )
  data.table::fwrite(out, out_path, sep = "\t", quote = FALSE, na = "NA",
                     compress = "gzip")
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

  # optional one-file spot check of the allele convention
  if (CHECK_N > 0 && !checked) {
    dt <- data.table::fread(mp, nrows = CHECK_N)
    cc <- intersect(c("SNP","CHR","BP","A1","A2","BETA","OR","FRQ"), names(dt))
    message("\n---- allele spot check (first ", CHECK_N, " rows of ",
            base, ") ----")
    message("   A1 = effect allele, FRQ = effect-allele frequency (recoded on input).")
    print(dt[, ..cc])
    message("---- confirm A1/FRQ match the cohort's effect_allele/effect_AF ----\n")
    checked <- TRUE
  }
}

message(sprintf("\nDONE.  %d munged, %d skipped, %d failed.  Output -> %s/",
                n_ok, n_skip, n_fail, OUT))
if (n_fail > 0) message("See per-file logs under ", log_root, "/ for failures.")
