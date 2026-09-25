#!/usr/bin/env Rscript
# ------------------------------------------------------------
# run_trecase.R
#
# Command line wrapper for asSeq::trecase.
#
# Usage:
#   Rscript run_trecase.R [options] <total.txt> <hap1.txt> <hap2.txt> \
#       <covariates.txt> <genotypes.vcf[.gz]> <tss.bed>
#
# Run with --help for a description of inputs and options.
# ------------------------------------------------------------

usage <- function() {
  cat('
Usage: run_trecase.R [options] TOTAL HAP1 HAP2 COVARIATES VCF TSS_BED

Run TReCASE eQTL mapping (asSeq::trecase) from the command line.

Positional arguments:
  TOTAL        Gene x sample matrix of total (gene-level) read counts. Tab
               delimited; first column holds gene IDs and the header row holds
               sample IDs. Used to build Y.
  HAP1         Gene x sample matrix of haplotype 1 allele-specific counts.
               Same format as TOTAL. Used to build Y1.
  HAP2         Gene x sample matrix of haplotype 2 allele-specific counts.
               Same format as TOTAL. Used to build Y2.
  COVARIATES   Covariate x sample matrix. Tab delimited; first column holds
               covariate names and the header row holds sample IDs. Used to
               build X. Do not include an intercept.
  VCF          VCF (optionally gzipped) with phased GT calls for the variants
               to test. Used to build Z, mChr and mPos. Genotypes are coded
               0|0 -> 0, 0|1 -> 1, 1|0 -> 3, 1|1 -> 4.
  TSS_BED      BED file of transcription start sites (chrom, start, end, gene
               ID). The 4th column must match gene IDs in TOTAL. Used to build
               eChr and ePos (ePos = start + 1).

Genes are restricted to those present in TOTAL, HAP1, HAP2 and TSS_BED, and
samples to those present in all matrices and the VCF.

trecase options:
  --output-tag=STR        Prefix for output files [required]. Writes
                          STR_eqtl.txt, STR_freq.txt, and STR_eqtl_annotated.txt
                          (eqtl results with gene and marker IDs).
  --p-cut=NUM             Only save associations with p-value < p-cut. Required
                          unless --permute is given; with --permute, the
                          nominal scan is only run if --p-cut is supplied.
  --offset=FILE           Per-sample offset: two-column file of sample ID and
                          offset value (optional header)
  --min-AS-reads=INT      [default 5]
  --min-AS-sample=INT     [default 5]
  --min-n-het=INT         [default 5]
  --local-only=BOOL       Test local eQTL only [default TRUE]
  --local-distance=INT    [default 200000]
  --converge=NUM          [default 5e-5]
  --convergeGLM=NUM       [default 1e-8]
  --scoreTestP=NUM        [default 0.05]
  --transTestP=NUM        [default 0.05]
  --trace=INT             [default 1]
  --maxit=INT             [default 100]

Permutation options (asSeq::trecaseP):
  --permute               Estimate gene-level permutation p-values for the
                          best local variant of each gene and write
                          STR_perm.txt. Requires --local-only=TRUE.
  --np-max=INT            Maximum number of permutations [default 5000]
  --np=LIST               Comma-separated, ascending numbers of permutations
                          at which to check for early stopping
                          [default 20,100,500,1000,2500]
  --aim-p=LIST            Comma-separated, descending p-value thresholds
                          matching --np. Permutation stops for a gene once its
                          permutation p-value is confidently above aim-p[i]
                          after np[i] permutations [default 0.5,0.2,0.1,0.05,0.02]
  --confidence-p=NUM      Binomial p-value used for the early stopping decision
                          [default 0.01]
  --seed=INT              Random seed for permutations

  STR_perm.txt has one row per gene. For each of the models trec, ase,
  trecase and trec_trecase (TReC p-value for trans-eQTL, TReCASE otherwise) it
  reports the best marker, its nominal p-value (pval_*), permutation p-value
  (perP_*, the fraction of permutations with a better minimum p-value, which
  can be 0), number of permutations (nuse_*), and a Benjamini-Hochberg q-value
  across genes (qval_*) computed from (k + 1) / (n + 1), where k is the number
  of permutations with a better minimum p-value and n is nuse_*.

Input handling options:
  --unphased=STR          How to treat unphased heterozygous calls (0/1):
                          "error" or "drop" (drop the variant) [default error]
  --drop-missing=BOOL     Drop variants with missing genotypes instead of
                          stopping with an error [default FALSE]
  --drop-low-variance=BOOL
                          Drop genes/variants whose variance across samples is
                          below --converge (trecase stops on these otherwise)
                          [default TRUE]
  --exclude-chroms=LIST   Comma-separated chromosomes to exclude, for both
                          genes and variants (e.g. X,Y,M). Names are matched
                          after the conversion described below, so X, chrX
                          and 23 are equivalent.
  -h, --help              Show this message and exit

Chromosome names are converted to integers: a leading "chr" is removed and
X, Y and M/MT are coded as 23, 24 and 25. Genes and variants on other contigs
are dropped. With --local-only=TRUE, variants on chromosomes with no genes to
test are dropped before genotypes are parsed, since they cannot be local to
any gene. Options may be given as --name=value or --name value; option
names are case-insensitive and "." or "_" may be used in place of "-".
')
}

# ------------------------------------------------------------
# argument parsing
# ------------------------------------------------------------

defaults <- list(
  "output-tag"        = NULL,
  "p-cut"             = NULL,
  "offset"            = NULL,
  "min-as-reads"      = 5,
  "min-as-sample"     = 5,
  "min-n-het"         = 5,
  "local-only"        = TRUE,
  "local-distance"    = 200000,
  "converge"          = 5e-5,
  "convergeglm"       = 1e-8,
  "scoretestp"        = 0.05,
  "transtestp"        = 0.05,
  "trace"             = 1,
  "maxit"             = 100,
  "unphased"          = "error",
  "drop-missing"      = FALSE,
  "drop-low-variance" = TRUE,
  "permute"           = FALSE,
  "np-max"            = 5000,
  "np"                = "20,100,500,1000,2500",
  "aim-p"             = "0.5,0.2,0.1,0.05,0.02",
  "confidence-p"      = 0.01,
  "seed"              = NULL,
  "exclude-chroms"    = NULL
)

# options that may be given without a value
flags <- c("permute")

die <- function(...) {
  message("Error: ", ...)
  quit(save = "no", status = 1)
}

parseBool <- function(x, name) {
  v <- toupper(x)
  if (v %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO")) return(FALSE)
  die(sprintf("--%s expects TRUE or FALSE, got '%s'", name, x))
}

parseNum <- function(x, name) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) die(sprintf("--%s expects a number, got '%s'", name, x))
  v
}

parseNumList <- function(x, name) {
  v <- suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
  if (length(v) == 0 || any(is.na(v))) {
    die(sprintf("--%s expects a comma-separated list of numbers, got '%s'", name, x))
  }
  v
}

parseArgs <- function(args) {
  opts <- defaults
  pos <- character(0)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0)
    }
    if (startsWith(a, "--")) {
      a <- substring(a, 3)
      if (grepl("=", a, fixed = TRUE)) {
        key <- sub("=.*$", "", a)
        val <- sub("^[^=]*=", "", a)
      } else if (tolower(gsub("[._]", "-", a)) %in% flags) {
        key <- a
        val <- "TRUE"
      } else {
        key <- a
        if (i == length(args)) die(sprintf("option --%s requires a value", key))
        i <- i + 1
        val <- args[i]
      }
      key <- tolower(gsub("[._]", "-", key))
      if (!key %in% names(defaults)) die(sprintf("unknown option --%s", key))
      opts[key] <- list(val)
    } else {
      pos <- c(pos, a)
    }
    i <- i + 1
  }

  if (length(pos) != 6) {
    usage()
    die(sprintf("expected 6 positional arguments, got %d", length(pos)))
  }
  if (is.null(opts[["output-tag"]])) die("--output-tag is required")
  opts[["permute"]] <- parseBool(as.character(opts[["permute"]]), "permute")
  if (is.null(opts[["p-cut"]]) && !opts[["permute"]]) {
    die("--p-cut is required unless --permute is given")
  }

  if (!is.null(opts[["p-cut"]])) opts[["p-cut"]] <- parseNum(opts[["p-cut"]], "p-cut")
  if (!is.null(opts[["seed"]])) opts[["seed"]] <- parseNum(opts[["seed"]], "seed")
  for (k in c("min-as-reads", "min-as-sample", "min-n-het",
              "local-distance", "converge", "convergeglm", "scoretestp",
              "transtestp", "trace", "maxit", "np-max", "confidence-p")) {
    opts[[k]] <- parseNum(opts[[k]], k)
  }
  for (k in c("local-only", "drop-missing", "drop-low-variance")) {
    opts[[k]] <- parseBool(as.character(opts[[k]]), k)
  }
  if (!opts[["unphased"]] %in% c("error", "drop")) {
    die("--unphased must be 'error' or 'drop'")
  }
  if (!is.null(opts[["exclude-chroms"]])) {
    ex <- strsplit(opts[["exclude-chroms"]], ",", fixed = TRUE)[[1]]
    exInt <- chrToInt(trimws(ex))
    if (length(ex) == 0 || any(is.na(exInt))) {
      die(sprintf("--exclude-chroms: unrecognized chromosome in '%s'", opts[["exclude-chroms"]]))
    }
    opts[["exclude-chroms"]] <- exInt
  }
  opts[["np"]] <- parseNumList(opts[["np"]], "np")
  opts[["aim-p"]] <- parseNumList(opts[["aim-p"]], "aim-p")
  if (opts[["permute"]]) {
    if (!opts[["local-only"]]) die("--permute requires --local-only=TRUE")
    if (length(opts[["np"]]) != length(opts[["aim-p"]])) {
      die("--np and --aim-p must have the same length")
    }
    if (any(diff(opts[["np"]]) <= 0)) die("--np must be strictly ascending")
    if (any(diff(opts[["aim-p"]]) >= 0)) die("--aim-p must be strictly descending")
    if (opts[["np-max"]] <= max(opts[["np"]])) die("--np-max must be larger than max(--np)")
  }

  names(pos) <- c("total", "hap1", "hap2", "covariates", "vcf", "bed")
  for (f in pos) if (!file.exists(f)) die(sprintf("file not found: %s", f))
  if (!is.null(opts[["offset"]]) && !file.exists(opts[["offset"]])) {
    die(sprintf("file not found: %s", opts[["offset"]]))
  }

  list(files = as.list(pos), opts = opts)
}

# ------------------------------------------------------------
# input readers
# ------------------------------------------------------------

# Read a feature x sample matrix; returns a numeric matrix with dimnames.
readMatrix <- function(file, what) {
  d <- read.delim(file, row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE, comment.char = "")
  m <- as.matrix(d)
  if (!is.numeric(m)) die(sprintf("non-numeric values found in %s", what))
  if (anyDuplicated(colnames(m))) die(sprintf("duplicated sample IDs in %s", what))
  message(sprintf("Read %s: %d rows x %d samples", what, nrow(m), ncol(m)))
  m
}

# Convert chromosome names to integer codes, NA for unrecognized contigs.
chrToInt <- function(chr) {
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- toupper(chr)
  chr[chr == "X"] <- "23"
  chr[chr == "Y"] <- "24"
  chr[chr %in% c("M", "MT")] <- "25"
  out <- suppressWarnings(as.integer(chr))
  out[!grepl("^[0-9]+$", chr)] <- NA
  out
}

readTSS <- function(file) {
  lines <- readLines(file)
  lines <- lines[!grepl("^(#|track|browser)", lines) & nzchar(lines)]
  f <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(f) < 4)) die("TSS BED file must have at least 4 columns")
  bed <- data.frame(chrom = vapply(f, `[`, "", 1),
                    start = as.numeric(vapply(f, `[`, "", 2)),
                    gene  = vapply(f, `[`, "", 4),
                    stringsAsFactors = FALSE)
  if (anyDuplicated(bed$gene)) {
    dup <- unique(bed$gene[duplicated(bed$gene)])
    message(sprintf("Warning: %d genes have multiple TSSs in the BED file; using the first entry for each",
                    length(dup)))
    bed <- bed[!duplicated(bed$gene), ]
  }
  bed$chr <- chrToInt(bed$chrom)
  bed$pos <- bed$start + 1
  message(sprintf("Read TSS BED: %d genes", nrow(bed)))
  bed
}

# Read a VCF and return the phased genotype matrix Z (samples x variants)
# along with mChr and mPos.
# keepChr / dropChr: integer chromosome codes to keep / drop, or NULL.
readVCF <- function(file, samples, unphased, dropMissing, keepChr = NULL,
                    dropChr = NULL) {
  con <- gzfile(file, "r")
  nMeta <- 0
  repeat {
    l <- readLines(con, n = 1)
    if (length(l) == 0) { close(con); die("no #CHROM header line found in VCF") }
    if (startsWith(l, "#CHROM")) break
    nMeta <- nMeta + 1
  }
  close(con)

  header <- strsplit(sub("^#", "", l), "\t", fixed = TRUE)[[1]]
  if (length(header) < 10 || header[9] != "FORMAT") {
    die("VCF has no FORMAT/sample columns")
  }
  vcf <- read.delim(gzfile(file), skip = nMeta + 1, header = FALSE,
                    col.names = header, check.names = FALSE,
                    colClasses = "character", comment.char = "",
                    quote = "", na.strings = character(0))
  message(sprintf("Read VCF: %d variants x %d samples", nrow(vcf), length(header) - 9))

  keepSamples <- intersect(samples, header[-(1:9)])

  # biallelic variants on recognized chromosomes only
  chr <- chrToInt(vcf$CHROM)
  ok <- !is.na(chr) & !grepl(",", vcf$ALT, fixed = TRUE)
  if (any(!ok)) {
    message(sprintf("Dropped %d multi-allelic or non-standard contig variants", sum(!ok)))
  }
  vcf <- vcf[ok, , drop = FALSE]
  chr <- chr[ok]

  if (!is.null(keepChr) || !is.null(dropChr)) {
    ok <- !chr %in% dropChr
    if (!is.null(keepChr)) ok <- ok & chr %in% keepChr
    if (any(!ok)) {
      message(sprintf("Dropped %d variants on chromosomes that are excluded or have no genes to test",
                      sum(!ok)))
    }
    vcf <- vcf[ok, , drop = FALSE]
    chr <- chr[ok]
  }

  # locate GT within FORMAT
  fmt <- strsplit(vcf$FORMAT, ":", fixed = TRUE)
  gtIdx <- vapply(fmt, function(x) match("GT", x), 1L)
  if (any(is.na(gtIdx))) die("some VCF records have no GT field")

  gt <- as.matrix(vcf[, keepSamples, drop = FALSE])
  if (all(gtIdx == 1L)) {
    gt[] <- sub(":.*$", "", gt)
  } else {
    for (j in seq_len(ncol(gt))) {
      gt[, j] <- mapply(function(s, k) strsplit(s, ":", fixed = TRUE)[[1]][k],
                        gt[, j], gtIdx)
    }
  }

  codes <- c("0|0" = 0, "0|1" = 1, "1|0" = 3, "1|1" = 4)
  Z <- matrix(codes[gt], nrow = nrow(gt), dimnames = dimnames(gt))

  unphasedHet <- gt %in% c("0/1", "1/0")
  dim(unphasedHet) <- dim(gt)
  badUnphased <- rowSums(unphasedHet) > 0
  if (any(badUnphased)) {
    if (unphased == "error") {
      die(sprintf("%d variants have unphased heterozygous genotypes; use phased genotypes or --unphased=drop",
                  sum(badUnphased)))
    }
    message(sprintf("Dropped %d variants with unphased heterozygous genotypes", sum(badUnphased)))
  }
  # homozygous unphased calls are unambiguous
  Z[gt == "0/0"] <- 0
  Z[gt == "1/1"] <- 4

  badMissing <- rowSums(is.na(Z)) > 0 & !badUnphased
  if (any(badMissing)) {
    if (!dropMissing) {
      die(sprintf("%d variants have missing or unrecognized genotypes; use --drop-missing=TRUE to drop them",
                  sum(badMissing)))
    }
    message(sprintf("Dropped %d variants with missing genotypes", sum(badMissing)))
  }

  keep <- !badUnphased & !badMissing
  ids <- ifelse(vcf$ID == ".", paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, sep = ":"), vcf$ID)
  Z <- t(Z[keep, , drop = FALSE])
  colnames(Z) <- ids[keep]

  list(Z = Z, mChr = chr[keep], mPos = as.numeric(vcf$POS[keep]),
       mID = ids[keep])
}

readOffset <- function(file) {
  d <- read.delim(file, header = FALSE, stringsAsFactors = FALSE,
                  colClasses = "character")
  if (ncol(d) < 2) die("offset file must have two columns: sample ID and offset")
  if (is.na(suppressWarnings(as.numeric(d[1, 2])))) d <- d[-1, , drop = FALSE]
  off <- as.numeric(d[, 2])
  if (any(is.na(off))) die("non-numeric values in offset file")
  names(off) <- d[, 1]
  off
}

# ------------------------------------------------------------
# analyses
# ------------------------------------------------------------

runNominal <- function(Y, Y1, Y2, X, Z, offset, eChr, ePos, mChr, mPos,
                       genes, mID, opts) {
  res <- trecase(Y, Y1, Y2, X, Z,
                 output.tag     = opts[["output-tag"]],
                 p.cut          = opts[["p-cut"]],
                 offset         = offset,
                 min.AS.reads   = opts[["min-as-reads"]],
                 min.AS.sample  = opts[["min-as-sample"]],
                 min.n.het      = opts[["min-n-het"]],
                 local.only     = opts[["local-only"]],
                 local.distance = opts[["local-distance"]],
                 eChr = eChr, ePos = ePos, mChr = mChr, mPos = mPos,
                 converge       = opts[["converge"]],
                 convergeGLM    = opts[["convergeglm"]],
                 scoreTestP     = opts[["scoretestp"]],
                 transTestP     = opts[["transtestp"]],
                 trace          = opts[["trace"]],
                 maxit          = opts[["maxit"]])

  # annotate results with gene and marker IDs
  eqtlFile <- sprintf("%s_eqtl.txt", opts[["output-tag"]])
  if (file.exists(eqtlFile)) {
    eqtl <- read.delim(eqtlFile, check.names = FALSE, colClasses = "character")
    annot <- data.frame(GeneID   = genes[as.integer(eqtl$GeneRowID)],
                        MarkerID = mID[as.integer(eqtl$MarkerRowID)],
                        eqtl, check.names = FALSE, stringsAsFactors = FALSE)
    write.table(annot, sprintf("%s_eqtl_annotated.txt", opts[["output-tag"]]),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  failed <- sum(res$yFailBaselineModel != 0)
  if (failed > 0) {
    message(sprintf("Baseline model failed for %d genes", failed))
  }
  if (res$succeed != 1) die("trecase did not complete successfully")
}

runPermutation <- function(Y, Y1, Y2, X, Z, offset, eChr, ePos, mChr, mPos,
                           genes, mID, opts) {
  if (!is.null(opts[["seed"]])) set.seed(opts[["seed"]])

  perm <- trecaseP(Y, Y1, Y2, X, Z,
                   offset         = offset,
                   min.AS.reads   = opts[["min-as-reads"]],
                   min.AS.sample  = opts[["min-as-sample"]],
                   min.n.het      = opts[["min-n-het"]],
                   local.only     = opts[["local-only"]],
                   local.distance = opts[["local-distance"]],
                   eChr = eChr, ePos = ePos, mChr = mChr, mPos = mPos,
                   converge       = opts[["converge"]],
                   convergeGLM    = opts[["convergeglm"]],
                   scoreTestP     = opts[["scoretestp"]],
                   transTestP     = opts[["transtestp"]],
                   np.max         = opts[["np-max"]],
                   np             = opts[["np"]],
                   aim.p          = opts[["aim-p"]],
                   confidence.p   = opts[["confidence-p"]],
                   trace          = opts[["trace"]],
                   maxit          = opts[["maxit"]])

  out <- data.frame(GeneID = genes[perm$geneID], stringsAsFactors = FALSE)
  for (type in c("trec", "ase", "trecase", "trec_trecase")) {
    perP <- perm[[paste0("perP_", type)]]
    nuse <- perm[[paste0("nuse_", type)]]
    # add-one estimate so that no gene gets a permutation p-value of 0
    k <- round(perP * nuse)
    out[[paste0("MarkerID_", type)]] <- mID[perm[[paste0("markerID_", type)]]]
    out[[paste0("pval_", type)]] <- perm[[paste0("pval_", type)]]
    out[[paste0("perP_", type)]] <- perP
    out[[paste0("nuse_", type)]] <- nuse
    out[[paste0("qval_", type)]] <- p.adjust((k + 1) / (nuse + 1), method = "BH")
  }

  # pval_* is 9 (or negative) when a gene could not be tested
  for (type in c("trec", "ase", "trecase", "trec_trecase")) {
    pv <- paste0("pval_", type)
    out[[pv]][is.na(out[[paste0("perP_", type)]])] <- NA
  }

  write.table(out, sprintf("%s_perm.txt", opts[["output-tag"]]),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message(sprintf("Wrote permutation results for %d genes", nrow(out)))
}

# ------------------------------------------------------------
# main
# ------------------------------------------------------------

main <- function() {
  a <- parseArgs(commandArgs(trailingOnly = TRUE))
  files <- a$files
  opts <- a$opts

  suppressPackageStartupMessages(library(asSeq))

  Ym  <- readMatrix(files$total, "total expression")
  Y1m <- readMatrix(files$hap1, "haplotype 1 expression")
  Y2m <- readMatrix(files$hap2, "haplotype 2 expression")
  Xm  <- readMatrix(files$covariates, "covariates")
  tss <- readTSS(files$bed)

  offset <- NULL
  if (!is.null(opts[["offset"]])) offset <- readOffset(opts[["offset"]])

  # samples shared by all inputs
  samples <- Reduce(intersect, list(colnames(Ym), colnames(Y1m),
                                    colnames(Y2m), colnames(Xm)))
  if (!is.null(offset)) samples <- intersect(samples, names(offset))

  # genes shared by all expression matrices and the TSS file
  tss <- tss[!is.na(tss$chr), ]
  exChr <- opts[["exclude-chroms"]]
  if (!is.null(exChr)) {
    ex <- tss$chr %in% exChr
    if (any(ex)) message(sprintf("Dropped %d genes on excluded chromosomes", sum(ex)))
    tss <- tss[!ex, ]
  }
  genes <- Reduce(intersect, list(rownames(Ym), rownames(Y1m),
                                  rownames(Y2m), tss$gene))
  if (length(genes) == 0) die("no genes in common across expression matrices and TSS BED")
  message(sprintf("Using %d genes present in all expression matrices and the TSS BED", length(genes)))

  # with local-only testing, variants can only be tested on chromosomes
  # that have genes
  keepChr <- NULL
  if (opts[["local-only"]]) keepChr <- unique(tss$chr[match(genes, tss$gene)])

  geno <- readVCF(files$vcf, samples, opts[["unphased"]], opts[["drop-missing"]],
                  keepChr, exChr)
  samples <- rownames(geno$Z)
  if (length(samples) == 0) die("no samples in common across all inputs")
  if (ncol(geno$Z) == 0) die("no variants left to test")
  message(sprintf("Using %d samples present in all inputs", length(samples)))

  Y  <- t(Ym[genes, samples, drop = FALSE])
  Y1 <- t(Y1m[genes, samples, drop = FALSE])
  Y2 <- t(Y2m[genes, samples, drop = FALSE])
  X  <- t(Xm[, samples, drop = FALSE])
  Z  <- geno$Z
  tss <- tss[match(genes, tss$gene), ]
  eChr <- as.numeric(tss$chr)
  ePos <- as.numeric(tss$pos)
  mChr <- as.numeric(geno$mChr)
  mPos <- geno$mPos
  mID  <- geno$mID
  if (!is.null(offset)) offset <- offset[samples]

  if (opts[["drop-low-variance"]]) {
    cv <- opts[["converge"]]
    lowY <- apply(Y, 2, var) < cv
    if (any(lowY)) {
      message(sprintf("Dropped %d genes with near-zero variance in total expression", sum(lowY)))
      Y <- Y[, !lowY, drop = FALSE]; Y1 <- Y1[, !lowY, drop = FALSE]
      Y2 <- Y2[, !lowY, drop = FALSE]
      eChr <- eChr[!lowY]; ePos <- ePos[!lowY]; genes <- genes[!lowY]
    }
    lowZ <- apply(Z, 2, var) < cv
    if (any(lowZ)) {
      message(sprintf("Dropped %d variants with near-zero genotype variance", sum(lowZ)))
      Z <- Z[, !lowZ, drop = FALSE]
      mChr <- mChr[!lowZ]; mPos <- mPos[!lowZ]; mID <- mID[!lowZ]
    }
    if (ncol(Y) == 0) die("no genes left after variance filtering")
    if (ncol(Z) == 0) die("no variants left after variance filtering")
  }
  message(sprintf("Testing %d genes and %d variants", ncol(Y), ncol(Z)))

  if (!is.null(opts[["p-cut"]])) {
    runNominal(Y, Y1, Y2, X, Z, offset, eChr, ePos, mChr, mPos, genes, mID, opts)
  }
  if (opts[["permute"]]) {
    runPermutation(Y, Y1, Y2, X, Z, offset, eChr, ePos, mChr, mPos, genes, mID, opts)
  }
  message("Done")
}

main()
