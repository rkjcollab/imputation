#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(
  description = "Rscript sex_check.R run by create_initial_input.sh.")

parser$add_argument(
  "-s", "--sexcheck", help="PLINK2 .sexcheck file path (required)", required=TRUE)
parser$add_argument(
  "-m", "--min-male-xf", help="Threshold for minimum male XF (required)", required=TRUE)
parser$add_argument(
  "-f", "--max-female-xf", help="Threshold for maximum female XF (required)", required=TRUE)
parser$add_argument(
  "-l", "--log", help="Create initial input log file (required)", required=TRUE)

args <- parser$parse_args()

# Plot sex check ---------------------------------------------------------------

sexcheck <- read.delim(args$sexcheck)
min_male_xf <- as.numeric(args$min_male_xf)
max_female_xf <- as.numeric(args$max_female_xf)

sexcheck$sex <- ifelse(
  sexcheck$F > min_male_xf, "Male",ifelse(sexcheck$F < max_female_xf, "Female",NA))

female_xf <- sexcheck$F[sexcheck$sex == "Female"]
male_xf   <- sexcheck$F[sexcheck$sex == "Male"]

path <- dirname(args$sexcheck)
png(paste0(path, "/sexcheck_plot.png"), width = 800, height = 600)
par(mar = c(5, 4, 4, 12))

breaks <- hist(sexcheck$F, breaks = 50, plot = FALSE)$breaks
female_h <- hist(female_xf, breaks = breaks, plot = FALSE)
male_h   <- hist(male_xf,   breaks = breaks, plot = FALSE)
ymax <- max(c(female_h$counts, male_h$counts)) * 1.2

hist(female_xf,
     breaks = breaks,
     col = "tomato",
     border = "black",
     xlab = "F",
     main = "Histogram of F by Sex",
     ylim = c(0, ymax))
hist(male_xf,
     breaks = breaks,
     col = "lightblue",
     border = "black",
     add = TRUE,
     ylim = c(0, ymax))
par(xpd = FALSE)
abline(v = max_female_xf, col = "tomato", lty = 2)
abline(v = min_male_xf, col = "lightblue", lty = 2)

par(xpd = TRUE)
legend("right",
       inset = c(-0.25, 0), 
       legend = c("Female", "Male"),
       fill = c("tomato", "lightblue"))
dev.off()

# Calculate sex check ----------------------------------------------------------

sexcheck_err <- sexcheck[sexcheck$STATUS != "OK", ]
sexcheck_err$case <- NA
sexcheck_err$case[is.na(sexcheck_err$SNPSEX)] <- "outside_f_thresh"
sexcheck_err$case[sexcheck_err$PEDSEX == 1 & sexcheck_err$SNPSEX == 2] <- "mislab_male"
sexcheck_err$case[sexcheck_err$PEDSEX == 2 & sexcheck_err$SNPSEX == 1] <- "mislab_female"

sexcheck_err_1 <- sexcheck_err[sexcheck_err$case == "mislab_male", ]
log_entry_1 <- data.frame(
  Category = "Sex-Specific",
  Description = "Remove samples mislabeled as male",
  Samples = paste0("(", nrow(sexcheck_err_1), ")"),
  SNPs = "()",
  stringsAsFactors = FALSE)

sexcheck_err_2 <- sexcheck_err[sexcheck_err$case == "mislab_female", ]
log_entry_2 <- data.frame(
  Category = "Sex-Specific",
  Description = "Remove samples mislabeled as female",
  Samples = paste0("(", nrow(sexcheck_err_2), ")"),
  SNPs = "()",
  stringsAsFactors = FALSE)

sexcheck_err_3 <- sexcheck_err[sexcheck_err$case == "outside_f_thresh", ]
log_entry_3 <- data.frame(
  Category = "Sex-Specific",
  Description = "Remove samples outside X-chr F thresholds",
  Samples = paste0("(", nrow(sexcheck_err_3), ")"),
  SNPs = "()",
  stringsAsFactors = FALSE)

# Write out --------------------------------------------------------------------

# Plot saved above when made

out_path <- paste0(gsub("\\.sexcheck", "", args$sexcheck), "_rm.txt")
write.table(
  sexcheck_err[, c("X.FID", "IID")],
  file = out_path,
  sep="\t",
  quote=F,
  row.names=F,
  col.names=F)

write.table(
  rbind(log_entry_1, log_entry_2, log_entry_3),
  file = args$log,
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F,
  append = T)
