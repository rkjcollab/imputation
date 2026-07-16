library(testthat)

make_afreq <- function(path, ids, alt_freqs) {
  afreq <- data.frame(
    `#CHROM` = rep(1, length(ids)),
    ID = ids,
    REF = "A",
    ALT = "G",
    `PROVISIONAL_REF?` = "Y",
    ALT_FREQS = alt_freqs,
    OBS_CT = 200
  )
  write.table(afreq, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

test_that("filt_mono writes only monomorphic SNP IDs and logs their count", {
  test_dir <- tempfile("filt_mono_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  afreq_file <- file.path(test_dir, "tmp_mono.afreq")
  make_afreq(
    afreq_file,
    ids = c("snp1", "snp2", "snp3", "snp4"),
    alt_freqs = c(0, 0.25, 0, 0.5)
  )

  log_file <- file.path(test_dir, "log.tsv")
  file.create(log_file)

  filt_mono(afreq_file, log_file)

  actual_ids <- readLines(file.path(test_dir, "tmp_mono_rm.txt"))
  expect_equal(actual_ids, c("snp1", "snp3"))

  log_lines <- readLines(log_file)
  expect_equal(
    log_lines,
    "Pre-Filtering\tRemove monomorphic SNPs\t()\t(2)"
  )
})

test_that("filt_mono handles case of no monomorphic SNPs", {
  test_dir <- tempfile("filt_mono_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  afreq_file <- file.path(test_dir, "tmp_mono.afreq")
  make_afreq(
    afreq_file,
    ids = c("snp1", "snp2"),
    alt_freqs = c(0.1, 0.5)
  )

  log_file <- file.path(test_dir, "log.tsv")
  file.create(log_file)

  filt_mono(afreq_file, log_file)

  actual_ids <- readLines(file.path(test_dir, "tmp_mono_rm.txt"))
  expect_equal(actual_ids, character(0))

  log_lines <- readLines(log_file)
  expect_equal(
    log_lines,
    "Pre-Filtering\tRemove monomorphic SNPs\t()\t(0)"
  )
})

test_that("filt_mono appends its log entry after existing log content", {
  test_dir <- tempfile("filt_mono_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  afreq_file <- file.path(test_dir, "tmp_mono.afreq")
  make_afreq(afreq_file, ids = "snp1", alt_freqs = 0)

  log_file <- file.path(test_dir, "log.tsv")
  writeLines(
    "Pre-Filtering\tRemove SNPs with missingness > 20%\t()\t(5)",
    log_file
  )

  filt_mono(afreq_file, log_file)

  log_lines <- readLines(log_file)
  expect_equal(
    log_lines,
    c(
      "Pre-Filtering\tRemove SNPs with missingness > 20%\t()\t(5)",
      "Pre-Filtering\tRemove monomorphic SNPs\t()\t(1)"
    )
  )
})

test_that("filt_mono throws an error if the .afreq file is empty", {
  test_dir <- tempfile("filt_mono_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  afreq_file <- file.path(test_dir, "tmp_mono.afreq")
  file.create(afreq_file)

  log_file <- file.path(test_dir, "log.tsv")
  file.create(log_file)

  expect_error(
    filt_mono(afreq_file, log_file),
    "PLINK2 .afreq file is empty"
  )
})

test_that("filt_mono throws an error if the ALT_FREQS column is missing", {
  test_dir <- tempfile("filt_mono_")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE))

  afreq_file <- file.path(test_dir, "tmp_mono.afreq")
  afreq <- data.frame(
    `#CHROM` = 1,
    ID = "snp1",
    REF = "A",
    ALT = "G",
    OBS_CT = 200
  )
  write.table(afreq, afreq_file, sep = "\t", quote = FALSE, row.names = FALSE)

  log_file <- file.path(test_dir, "log.tsv")
  file.create(log_file)

  expect_error(
    filt_mono(afreq_file, log_file),
    "ALT_FREQS column is missing from .afreq file"
  )
})
