# ================================================================
# Phase Screen Score histogram + stats exporter
# - Reads 1.csv ... 13.csv (tab-delimited) with columns: Title, PhaseScreenScore
# - Creates histogram (unchanged aesthetics)
# - Saves plot as Fig2B.png in the same directory as this script
# - Exports summary statistics to Stats.xlsx (same directory)
# ================================================================

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(writexl)
  library(e1071)   # for skewness/kurtosis
})

# -----------------------------
# 1) Resolve "script directory"
#    Works in Source, Rscript, or Console (falls back to getwd()).
# -----------------------------
script_dir <- (function() {
  p1 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
  if (!is.null(p1) && nzchar(p1)) return(dirname(p1))
  p2 <- tryCatch(normalizePath(sub("^--file=", "",
                                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])),
                 error = function(e) NULL)
  if (!is.null(p2) && nzchar(p2)) return(dirname(p2))
  p3 <- tryCatch(normalizePath(getOption("knit.rmd.original", "")), error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3)) return(dirname(p3))
  getwd()
})()

# -----------------------------
# 2) Read and combine data (1.csv ... 13.csv)
#    - Uses tab delimiter
#    - Selects Title and PhaseScreenScore
#    - Drops the first row (slice(2:n())) as in your original code
# -----------------------------
read_one <- function(i) {
  fpath <- file.path(script_dir, sprintf("%d.csv", i))
  df <- read_delim(
    file = fpath,
    delim = "\t",
    show_col_types = FALSE,
    col_select = c("Title", "PhaseScreenScore")
  )
  if (nrow(df) >= 2) df <- dplyr::slice(df, 2:n())
  df
}

files_idx <- 1:13
megadata  <- dplyr::bind_rows(lapply(files_idx, read_one))

# Optional cleaning (kept separate, your plot uses 'megadata' directly)
megadata_clean <- megadata[is.finite(megadata$PhaseScreenScore), ]

# -----------------------------
# 3) Plot (keep EXACT aesthetics as provided)
# -----------------------------
cutoff <- 0.797

hypo_db <- ggplot(megadata, aes(x = PhaseScreenScore)) +
  geom_histogram(
    aes(fill = after_stat(x > cutoff)),
    position = "dodge",
    colour = "black",
    alpha = 1,
    bins = 50
  ) +
  scale_fill_manual(values = c("FALSE" = "#CAE1FF", "TRUE" = "#8B3A3A")) +
  labs(x = "Phase Screen Score", y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 1),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title   = element_text(color = "black", face = "bold", size = 22),
    axis.text.x  = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y  = element_text(size = 16, color = "black", face = "bold"),
    legend.position = "none",
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    strip.text = element_blank()
  ) +
  # geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1) +
  scale_y_continuous(labels = comma, breaks = seq(0, 200000, 50000)) +
  coord_flip()

# Preview (optional)
print(hypo_db)

# Save plot as Fig2B.png in the same directory as this script
ggsave(
  filename = file.path(script_dir, "Fig2B.png"),
  plot = hypo_db,
  width = 7,
  height = 5,
  dpi = 600,
  units = "in",
  bg = "white"
)

# -----------------------------
# 4) Summary statistics and export to Excel
# -----------------------------
results <- megadata %>%
  summarise(
    Min       = round(min(PhaseScreenScore, na.rm = TRUE), 3),
    Q1        = round(quantile(PhaseScreenScore, 0.25, na.rm = TRUE), 3),
    Median    = round(median(PhaseScreenScore, na.rm = TRUE), 3),
    Q3        = round(quantile(PhaseScreenScore, 0.75, na.rm = TRUE), 3),
    Max       = round(max(PhaseScreenScore, na.rm = TRUE), 3),
    Mean      = round(mean(PhaseScreenScore, na.rm = TRUE), 3),
    SD        = round(sd(PhaseScreenScore, na.rm = TRUE), 3),
    Variance  = round(var(PhaseScreenScore, na.rm = TRUE), 3),
    Range     = round(max(PhaseScreenScore, na.rm = TRUE) - min(PhaseScreenScore, na.rm = TRUE), 3),
    IQR       = round(IQR(PhaseScreenScore, na.rm = TRUE), 3),
    Skewness  = round(e1071::skewness(PhaseScreenScore, na.rm = TRUE), 3),
    Kurtosis  = round(e1071::kurtosis(PhaseScreenScore, na.rm = TRUE), 3)
  )

# Rows strictly greater than Q3 (as in your original logic)
rows_above_Q3 <- megadata %>%
  filter(PhaseScreenScore > results$Q3)

n_rows_above_Q3 <- nrow(rows_above_Q3)

# Prepare an additional small table for the Q3 count
q3_info <- tibble::tibble(
  Q3 = results$Q3,
  Count_above_Q3 = n_rows_above_Q3
)

# Write to Excel (Stats.xlsx) in the same directory
write_xlsx(
  x = list(PhaseScreenScore_stats = results,
           Above_Q3 = q3_info),
  path = file.path(script_dir, "Stats.xlsx")
)

# -----------------------------
# 5) Done
# -----------------------------
message("Saved plot:  ", file.path(script_dir, "Fig2B.png"))
message("Saved stats: ", file.path(script_dir, "Stats.xlsx"))
message("Count of rows above Q3: ", n_rows_above_Q3)
