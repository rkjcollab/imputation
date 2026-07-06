#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

print("First trailing arg should be pre_qc directory path")
print("Second trailing arg should be post_qc directory path")
print("This _updated version fixes both strand flips and")
print("strand flip and allele switch. But, in the case of")
print("duplicate variants, or duplicate chr:pos locations")
print("where one has an allele mismatch, both variants are")
print("removed due to inability to identify which of the")
print("pair is correct.")

# TODO: when revisit and rewrite, remove imp_server arg

get_strand_flip_snp_names <- function(pre_qc_dir, post_qc_dir, imp_server) {
  # Load snp-excluded.txt
  snp_excl <- read.delim(paste0(pre_qc_dir, "/snps-excluded.txt"),
    stringsAsFactors = FALSE
  )

  # Add column with chr:pos for all variants
  # TODO: revisit and simplify now that both match
  snp_excl$chr.pos <- paste0(snp_excl$CHROM, ":", snp_excl$POS)
  snp_excl$FilterType <- snp_excl$INFO # copy to same col name as tm
  snp_excl$Info <- snp_excl$INFO # copy to same col name as tm
  snp_excl$ref <- snp_excl$REF # copy to same col name as tm

  # Load bim so can get actual varID names
  bim <- read.table(paste0(pre_qc_dir, "/pre_qc.bim"),
    stringsAsFactors = FALSE
  )
  bim$chr.pos <- paste0(bim$V1, ":", bim$V4)

  # Make sure presence/absence of chr matches bim file
  if (any(grepl("^chr", bim$V1))) {
    snp_excl$chr.pos <- ifelse(
      any(grepl("^chr", snp_excl$chr.pos)),
      snp_excl$chr.pos,
      paste0("chr", snp_excl$chr.pos)
    )
  } else {
    snp_excl$chr.pos <- gsub("^chr", "", snp_excl$chr.pos)
  }

  # Add varID names to snp_excl by chr:pos
  merged <- merge(snp_excl, bim, by = c("chr.pos"))

  # Get strand flips & strand flip and allele switch
  snp_frame_flip_as <- merged[grep("Strand flip", merged$FilterType), ]

  # Remove strand flip & strand flip and allele switch from from snp_excl
  merged_to_excl <-
    merged[!merged$FilterType %in% snp_frame_flip_as$FilterType, ]

  # Get allele switch only snps
  snp_frame_as <-
    merged_to_excl[grep("Allele switch", merged_to_excl$FilterType), ]

  ### Handle strand flips & strand flip and allele switches
  # For both strand flip & strand flip and allele switch, flip strand
  # Only do if strand flips were found, otherwise create empty list
  if (nrow(snp_frame_flip_as) > 0) {
    # Get varIDs from .bim by chr:pos for strand flip & strand flip and
    # allele switch
    snps <- bim$V2[bim$chr.pos %in% snp_frame_flip_as$chr.pos]

    # Get ref alleles for strand flip and allele switch
    snp_frame_flip_as_both <-
      snp_frame_flip_as[
        grep("Strand flip and Allele switch", snp_frame_flip_as$FilterType),
      ]
    snp_frame_flip_as_both$ref <- gsub(
      "/[[:alpha:]]", "",
      gsub(".*:", "", snp_frame_flip_as_both$Info)
    )

    # Get varID name with reference a2 allele for PLINK --a2-allele
    snp_frame_flip_as_both <- snp_frame_flip_as_both[, c("V2", "ref")]
  } else {
    snps <- character()
    snp_frame_flip_as_both <- character()
    print("No strand flip SNPs were found in the excluded SNPs.")
  }

  # Write out SNPs for input into PLINK --flip, will flip both strand
  # flip & strand flip and allele switch
  write.table(snps, paste0(post_qc_dir, "/tmp_flip.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  # Write out SNPs for input into PLINK --a2-allele, will allele switch
  # strand flip and allele switch alleles after --flip
  write.table(
    snp_frame_flip_as_both, paste0(post_qc_dir, "/tmp_a2-allele.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  ### Handle allele switches only
  # Get ref alleles for strand flip and allele switch
  snp_frame_as$ref <- gsub(
    "/[[:alpha:]]", "",
    gsub(".*:", "", snp_frame_as$Info)
  )

  # Get varID name with reference a2 allele for PLINK --a2-allele
  snp_frame_as <- snp_frame_as[, c("V2", "ref")]

  # Write out SNPs for input into PLINK --a2-allele, will allele switch
  # alleles switches only
  write.table(
    snp_frame_as, paste0(post_qc_dir, "/tmp_a2-allele_switch_only.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
}

get_strand_flip_snp_names(args[1], args[2], args[3])
