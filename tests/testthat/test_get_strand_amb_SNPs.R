library(testthat)

test_that("get_strand_amb_SNPs writes only strand ambiguous SNP IDs", {
  test_dir <- tempfile("strand_amb_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  # Make BIM file with strand ambiguous SNPs
  bim_file <- file.path(test_dir, ".bim")
  bim <- data.frame(
    chr = rep(1, 6),
    id = c("at_snp", "ta_snp", "cg_snp", "gc_snp", "ac_snp", "ga_snp"),
    cm = 0,
    pos = seq(100, 600, 100),
    a1 = c("A", "T", "C", "G", "A", "G"),
    a2 = c("T", "A", "G", "C", "C", "A")
  )
  write.table(bim, bim_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )

  get_strand_amb_SNPs(bim_file)

  actual <- readLines(file.path(test_dir, "tmp_strand_remove_snps.txt"))
  expected <- c("at_snp", "ta_snp", "cg_snp", "gc_snp")
  expect_equal(actual, expected)
})

test_that(
  "get_strand_amb_SNPs writes empty file when no strand ambiguous SNPs",
  {
    test_dir <- tempfile("no_strand_amb_")
    dir.create(test_dir, recursive = TRUE)
    on.exit(unlink(test_dir, recursive = TRUE))

    # Make BIM file with no strand ambiguous SNPs
    bim_file <- file.path(test_dir, ".bim")
    bim <- data.frame(
      chr = rep(1, 4),
      id = c("ac_snp", "ga_snp", "ca_snp", "ag_snp"),
      cm = 0,
      pos = seq(100, 400, 100),
      a1 = c("A", "G", "C", "A"),
      a2 = c("C", "A", "A", "G")
    )
    write.table(bim, bim_file,
      sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE
    )

    get_strand_amb_SNPs(bim_file)

    actual <- readLines(file.path(test_dir, "tmp_strand_remove_snps.txt"))
    expected <- character(0)
    expect_equal(actual, expected)
  }
)

test_that("get_strand_amb_SNPs throws an error when BIM file is empty", {
  test_dir <- tempfile("empty_bim_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  # Make empty BIM file
  bim_file <- file.path(test_dir, ".bim")
  write.table(data.frame(), bim_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )

  expect_error(get_strand_amb_SNPs(bim_file), "BIM file is empty")
})
