#!/usr/bin/env Rscript

# Script is called by create_initial_input. Requires PLINK2 .afreq file and
# path to log file. Removes monomorphic SNPs.

library(argparse)

filt_mono <- function(afreq_file, log_file) {
  # Throw error if .afreq file is empty
  if (file.size(afreq_file) == 0) {
    stop("PLINK2 .afreq file is empty")
  }

  afreq <- read.delim(afreq_file)

  # Throw error if ALT_FREQS column is missing
  if (!"ALT_FREQS" %in% colnames(afreq)) {
    stop("ALT_FREQS column is missing from .afreq file")
  }

  afreq_rm <- afreq[afreq$ALT_FREQS == 0, ]
  afreq_rm_ct <- nrow(afreq_rm)

  log_entry <- data.frame(
    Category = "Pre-Filtering",
    Description = "Remove monomorphic SNPs",
    Samples = "()",
    SNPs = paste0("(", afreq_rm_ct, ")"),
    stringsAsFactors = FALSE
  )

  path <- dirname(afreq_file)
  write.table(
    afreq_rm$ID,
    file = paste0(path, "/tmp_mono_rm.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  write.table(
    log_entry,
    file = log_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
}

main <- function() {
  parser <- ArgumentParser(
    description = "Rscript filt_mono.R run by create_initial_input.sh."
  )

  parser$add_argument(
    "-a", "--afreq",
    help = "PLINK2 .afreq file path (required)", required = TRUE
  )
  parser$add_argument(
    "-l", "--log",
    help = "Create initial input log file (required)", required = TRUE
  )

  args <- parser$parse_args()
  filt_mono(args$afreq, args$log)
}

if (sys.nframe() == 0) {
  main()
}
