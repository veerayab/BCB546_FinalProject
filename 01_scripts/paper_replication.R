############################
# Reproduce Figure 1A, Table 1, Table S3
# Paper: Himes et al. 2014
# GEO: GSE52778
#
# This script reproduces only the RNA-seq differential
# expression outputs that are available from GEO.
# It does NOT reproduce the wet-lab validation figures.
############################

############################
# 0) Install and load packages
############################

pkgs <- c("dplyr", "ggplot2", "tibble")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(ggplot2)
library(tibble)

############################
# 1) Create directories
############################

dir.create("00_data", showWarnings = FALSE)
dir.create("01_scripts", showWarnings = FALSE)

dir.create("02_output_files", showWarnings = FALSE)
dir.create("02_output_files/results", showWarnings = FALSE)
dir.create("02_output_files/figures", showWarnings = FALSE)
dir.create("02_output_files/tables", showWarnings = FALSE)

############################
# 2) Download GEO processed files
############################

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/suppl"

fpkm_file <- file.path("00_data", "GSE52778_All_Sample_FPKM_Matrix.txt.gz")
diff_file <- file.path("00_data", "GSE52778_Dex_vs_Untreated_gene_exp.diff.gz")

if (!file.exists(fpkm_file)) {
  download.file(
    url = paste0(base_url, "/GSE52778_All_Sample_FPKM_Matrix.txt.gz"),
    destfile = fpkm_file,
    mode = "wb"
  )
}

if (!file.exists(diff_file)) {
  download.file(
    url = paste0(base_url, "/GSE52778_Dex_vs_Untreated_gene_exp.diff.gz"),
    destfile = diff_file,
    mode = "wb"
  )
}

############################
# 3) Read files
############################

fpkm <- read.delim(
  gzfile(fpkm_file),
  check.names = TRUE,
  stringsAsFactors = FALSE
)

diff <- read.delim(
  gzfile(diff_file),
  check.names = TRUE,
  stringsAsFactors = FALSE
)

############################
# 4) Inspect column names
############################

cat("FPKM file columns:\n")
print(names(fpkm))

cat("\nDiff file columns:\n")
print(names(diff))

############################
# 5) Auto-detect important columns
############################

gene_col <- grep("^gene$|gene_short_name", names(diff), ignore.case = TRUE, value = TRUE)[1]
locus_col <- grep("^locus$", names(diff), ignore.case = TRUE, value = TRUE)[1]
testid_col <- grep("^test_id$", names(diff), ignore.case = TRUE, value = TRUE)[1]
sample1_col <- grep("^sample_1$", names(diff), ignore.case = TRUE, value = TRUE)[1]
sample2_col <- grep("^sample_2$", names(diff), ignore.case = TRUE, value = TRUE)[1]
value1_col <- grep("^value_1$", names(diff), ignore.case = TRUE, value = TRUE)[1]
value2_col <- grep("^value_2$", names(diff), ignore.case = TRUE, value = TRUE)[1]
teststat_col <- grep("^test_stat$", names(diff), ignore.case = TRUE, value = TRUE)[1]
sig_col <- grep("^significant$", names(diff), ignore.case = TRUE, value = TRUE)[1]
status_col <- grep("^status$", names(diff), ignore.case = TRUE, value = TRUE)[1]
p_col <- grep("^p[._]?value$", names(diff), ignore.case = TRUE, value = TRUE)[1]
q_col <- grep("^q[._]?value$", names(diff), ignore.case = TRUE, value = TRUE)[1]
fc_col <- grep("log2.*fold.*change", names(diff), ignore.case = TRUE, value = TRUE)[1]


############################
# 6) Clean the differential expression table
############################
# Why this is needed:
# The raw diff file may contain absurd placeholder values
# such as -1.79769e+308 for log2 fold change.
# Those break the volcano plot by stretching the x-axis.

clean_df <- diff %>%
  mutate(
    log2fc = suppressWarnings(as.numeric(.data[[fc_col]])),
    pval   = suppressWarnings(as.numeric(.data[[p_col]])),
    qval   = suppressWarnings(as.numeric(.data[[q_col]])),
    value1_num = suppressWarnings(as.numeric(.data[[value1_col]])),
    value2_num = suppressWarnings(as.numeric(.data[[value2_col]]))
  )

# Keep only rows with status OK, if the status column exists
if (!is.na(status_col) && !is.null(status_col) && nzchar(status_col)) {
  clean_df <- clean_df %>% filter(.data[[status_col]] == "OK")
}

# Remove rows with non-finite or impossible values
clean_df <- clean_df %>%
  filter(is.finite(log2fc), is.finite(pval), is.finite(qval)) %>%
  filter(pval > 0, qval >= 0) %>%
  filter(abs(log2fc) < 50) %>%   # removes absurd placeholder values
  mutate(sig = qval < 0.05)

############################
# 7) Sanity checks
############################

cat("\nSummary of cleaned log2 fold change:\n")
print(summary(clean_df$log2fc))

cat("\nRange of cleaned log2 fold change:\n")
print(range(clean_df$log2fc, na.rm = TRUE))

cat("\nNumber of significant genes with q < 0.05:\n")
print(sum(clean_df$sig, na.rm = TRUE))

############################
# 8) Figure 1A: Volcano plot
############################

p_volcano <- ggplot(clean_df, aes(x = log2fc, y = -log10(pval))) +
  geom_point(aes(color = sig), size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue")) +
  coord_cartesian(xlim = c(-15, 15), ylim = c(0, 16)) +
  labs(
    x = expression(Log[2] * "[Fold Change]"),
    y = expression(-Log[10] * "[P-value]"),
    title = "Figure 1A replication: Dex vs Untreated"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

print(p_volcano)

ggsave(
  filename = "02_output_files/figures/Figure1A_volcano_replication_fixed.png",
  plot = p_volcano,
  width = 5.5,
  height = 5,
  dpi = 300
)

############################
# 9) Table 1 replication
############################
# Paper Table 1 = top genes with Q-value < 1E-10

table1_rep <- clean_df %>%
  filter(qval < 1e-10) %>%
  transmute(
    Gene = .data[[gene_col]],
    locus = .data[[locus_col]],
    Mean_FPKM_Control = value1_num,
    Mean_FPKM_Dex = value2_num,
    log2_Fold_Change = log2fc,
    P_value = pval,
    Q_value = qval
  ) %>%
  arrange(Q_value, P_value, desc(abs(log2_Fold_Change)))

cat("\nTable 1 replication:\n")
print(table1_rep)
tibble::as_tibble(table1_rep) |> print(n = Inf)
View(table1_rep)

write.csv(
  table1_rep,
  file = "02_output files/tables/Table1_replication_fixed.csv",
  row.names = FALSE
)

############################
# 10) Table S3 replication
############################
# Paper Table S3 = all significant genes with q < 0.05

tableS3_rep <- clean_df %>%
  filter(qval < 0.05) %>%
  transmute(
    test_id = .data[[testid_col]],
    gene = .data[[gene_col]],
    locus = .data[[locus_col]],
    sample_1 = .data[[sample1_col]],
    sample_2 = .data[[sample2_col]],
    value_1 = value1_num,
    value_2 = value2_num,
    log2_fold_change = log2fc,
    test_stat = .data[[teststat_col]],
    p_value = pval,
    q_value = qval,
    significant = .data[[sig_col]]
  ) %>%
  arrange(q_value, p_value)

cat("\nNumber of rows in Table S3 replication:\n")
print(nrow(tableS3_rep))
view(tableS3_rep)
write.csv(
  tableS3_rep,
  file = "02_output files/tables/TableS3_replication_fixed.csv",
  row.names = FALSE
)

############################
# 11) Session info
############################

sessionInfo()
