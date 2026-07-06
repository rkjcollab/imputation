library(testthat)

test_that("get_strand_amb_SNPs writes only strand ambiguous SNP IDs", {
  test_dir <- tempfile("strand_amb_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  bim_file <- file.path(test_dir, "input.bim")
  bim <- data.frame(
    V1 = rep(1, 6),
    V2 = c("at_snp", "ta_snp", "cg_snp", "gc_snp", "ac_snp", "ga_snp"),
    V3 = 0,
    V4 = seq(100, 600, 100),
    V5 = c("A", "T", "C", "G", "A", "G"),
    V6 = c("T", "A", "G", "C", "C", "A")
  )

  write.table(bim, bim_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  get_strand_amb_SNPs(bim_file)

  actual <- readLines(file.path(test_dir, "tmp_strand_remove_snps.txt"))
  expected <- c("at_snp", "ta_snp", "cg_snp", "gc_snp")

  expect_equal(actual, expected)
})
