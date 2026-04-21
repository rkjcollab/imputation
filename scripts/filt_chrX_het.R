#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(
  description = "Rscript filt_chrX_het.R run by create_initial_input.sh.")

# parser$add_argument(
#   "-h", "--hetx", help="PLINK .hardy.x file path (required)", required=TRUE)
parser$add_argument(
  "-i", "--hh", help="PLINK1.9 .hh file path (required)", required=TRUE)
parser$add_argument(
  "-m", "--n-male", help="Number males in dataset (required)", required=TRUE)
parser$add_argument(
  "-p", "--perc-het", help="Percent (0-100) heterozygosity threshold (required)", required=TRUE)
parser$add_argument(
  "-l", "--log", help="Create initial input log file (required)", required=TRUE)

args <- parser$parse_args()

# chrX Heterozygosity ----------------------------------------------------------

# Use PLINK1.9's .hh file:
  # Produced automatically when the input data contains heterozygous calls where
  # they shouldn't be possible (haploid chromosomes, male X/Y), or there are
  # nonmissing calls for nonmales on the Y chromosome. A text file with one line
  # per error (sorted primarily by variant ID, secondarily by sample ID) with
  # the following three fields: FID, IID, variant ID

hh <- read.delim(args$hh, col.names = c("FID", "IID", "ID"))
n_male <- as.numeric(args$n_male)
perc_het <- as.numeric(args$perc_het)

id_count <- as.data.frame(table(hh$ID), stringsAsFactors = FALSE)
names(id_count) <- c("ID", "ID_count")
id_count$ID_count <- as.integer(id_count$ID_count)

hh_summ <- data.frame(
  ID = id_count$ID,
  ID_count = id_count$ID_count,
  het_perc = id_count$ID_count / n_male,
  stringsAsFactors = FALSE
)

hh_summ_rm <- hh_summ[hh_summ$het_perc > (perc_het / 100), , drop = FALSE]
hh_summ_rm_ct = length(unique(hh_summ_rm$ID))

log_entry <- data.frame(
  Category = "Sex-Specific",
  Description = paste0("Remove X-chr SNPs with heterozygosity >", perc_het, "% in males"),
  Samples = "()",
  SNPs = paste0("(", hh_summ_rm_ct, ")"),
  stringsAsFactors = FALSE)

# Write out --------------------------------------------------------------------

out_path <- paste0(gsub("\\.hh", "", args$hh), ".txt")
write.table(
  hh_summ_rm$ID,
  file = out_path,
  sep="\t",
  quote=F,
  row.names=F,
  col.names=F)

write.table(
  log_entry,
  file = args$log,
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F,
  append = T)
