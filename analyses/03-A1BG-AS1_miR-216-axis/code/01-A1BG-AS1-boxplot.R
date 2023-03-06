# Differential Expression Analysis of GNG C1 vs C2 using DESeq2
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html)

## LOAD LIBRARIES
library(dplyr)
library(ggpubr)

## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "03-A1BG-AS1_miR-216-axis")

# Set output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Set input directory
data_dir <- file.path(root_dir, "data")

# Declare input file paths
data_file <- file.path(data_dir, "pbta_raw_counts_mrna.rds")
metadata_file <- file.path(data_dir, "pbta_meta_mrna.rds")
# Declare output file paths
boxplot_pdf <- file.path(plots_dir, "a1bg-as1_boxplot.pdf")

#######
# Import metadata and data 
metadata <- readr::read_rds(metadata_file)
expression_df <- readr::read_rds(data_file)

#######
## PREPROCESS THE DATA
df <- expression_df %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample_id")  %>%
  inner_join(metadata, by="sample_id") %>%
  filter(!short_histology %in% c('GNT', 'Schwannoma', 'Teratoma'))


## Make box plot with stats
bp <- ggboxplot(df, 
  # Specify x values
  x = "short_histology",
  # Specify y values
  y = "`A1BG-AS1`",
  # Color in the box plot
  fill = "short_histology",
  # Specify color palette
  palette = "jco", 
  # Add x-axis label
  xlab="Brain Tumor Type",
  # Add y-axis label
  ylab="A1BG-AS1 Gene Expression",
  # Add points
  add = "jitter") +
  # Add p-value
  stat_compare_means()

# View boxplot
print(bp)

# Export plot
ggsave(boxplot_pdf, bp, width = 10, height = 6)

dev.off()

