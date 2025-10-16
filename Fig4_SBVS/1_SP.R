# ================================================================
# SP docking: geometric-mean summaries + figures (unchanged style)
# - Reads: combined_data_SP.rds, controls_SP.csv
# - Outputs: Histogram.png, PerFrame.png, Stats_SP.xlsx
# - All files saved in the same directory as this script
# - Plot aesthetics are preserved exactly as in your original code
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(e1071)   # skewness/kurtosis
  library(writexl) # write Excel with multiple sheets
})

# -----------------------------
# 1) Resolve "script directory"
#    Works in Source, Rscript, or Console (fallback to getwd()).
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
# 2) Load inputs
# -----------------------------
combined_data_SP <- readRDS(file.path(script_dir, "combined_data_SP.rds")) %>%
  dplyr::rename(Ligand = Title, Energy = docking_score)

controls_data_SP <- fread(file.path(script_dir, "controls_SP.csv")) %>%
  dplyr::rename(Ligand = Title, Energy = docking_score)

# Optional full dataset (not used downstream but kept for completeness)
# data_SP <- rbind(controls_data_SP, combined_data_SP)

# -----------------------------
# 3) Geometric mean (VS molecules and per frame)
#    NOTE: uses the same signed geometric-mean form as your original code
# -----------------------------
geo_mean_SP <- combined_data_SP %>%
  group_by(Ligand) %>%
  summarise(Geometric_mean = sign(mean(Energy)) * exp(mean(log(abs(Energy)))), .groups = "drop")

geo_mean_per_frame_SP <- combined_data_SP %>%
  group_by(Ligand, Frame) %>%
  summarise(Geometric_mean = sign(mean(Energy)) * exp(mean(log(abs(Energy)))), .groups = "drop")

# -----------------------------
# 4) Geometric mean for controls (overall and per frame)
# -----------------------------
geo_mean_control_SP <- controls_data_SP %>%
  group_by(Ligand) %>%
  summarise(Geometric_mean = sign(mean(Energy)) * exp(mean(log(abs(Energy)))), .groups = "drop")

geo_mean_control_per_frame_SP <- controls_data_SP %>%
  group_by(Ligand, Frame) %>%
  summarise(Geometric_mean = sign(mean(Energy)) * exp(mean(log(abs(Energy)))), .groups = "drop")

# -----------------------------
# 5) Summary stats (VS overall and per frame)
# -----------------------------
stat_vs_SP <- geo_mean_SP %>%
  summarise(
    Min      = round(min(Geometric_mean, na.rm = TRUE), 3),
    Q1       = round(quantile(Geometric_mean, 0.25, na.rm = TRUE), 3),
    Median   = round(median(Geometric_mean, na.rm = TRUE), 3),
    Q3       = round(quantile(Geometric_mean, 0.75, na.rm = TRUE), 3),
    Max      = round(max(Geometric_mean, na.rm = TRUE), 3),
    Mean     = round(mean(Geometric_mean, na.rm = TRUE), 3),
    SD       = round(sd(Geometric_mean, na.rm = TRUE), 3),
    Variance = round(var(Geometric_mean, na.rm = TRUE), 3),
    Range    = round(max(Geometric_mean, na.rm = TRUE) - min(Geometric_mean, na.rm = TRUE), 3),
    IQR      = round(IQR(Geometric_mean, na.rm = TRUE), 3),
    Skewness = round(e1071::skewness(Geometric_mean, na.rm = TRUE), 3),
    Kurtosis = round(e1071::kurtosis(Geometric_mean, na.rm = TRUE), 3)
  )

stat_vs_per_frame_SP <- geo_mean_per_frame_SP %>%
  group_by(Frame) %>%
  summarise(
    Min      = round(min(Geometric_mean, na.rm = TRUE), 3),
    Q1       = round(quantile(Geometric_mean, 0.25, na.rm = TRUE), 3),
    Median   = round(median(Geometric_mean, na.rm = TRUE), 3),
    Q3       = round(quantile(Geometric_mean, 0.75, na.rm = TRUE), 3),
    Max      = round(max(Geometric_mean, na.rm = TRUE), 3),
    Mean     = round(mean(Geometric_mean, na.rm = TRUE), 3),
    SD       = round(sd(Geometric_mean, na.rm = TRUE), 3),
    Variance = round(var(Geometric_mean, na.rm = TRUE), 3),
    Range    = round(max(Geometric_mean, na.rm = TRUE) - min(Geometric_mean, na.rm = TRUE), 3),
    IQR      = round(IQR(Geometric_mean, na.rm = TRUE), 3),
    Skewness = round(e1071::skewness(Geometric_mean, na.rm = TRUE), 3),
    Kurtosis = round(e1071::kurtosis(Geometric_mean, na.rm = TRUE), 3),
    .groups  = "drop"
  )

# -----------------------------
# 6) Summary stats (controls overall and per frame)
# -----------------------------
stat_controls_SP <- controls_data_SP %>%
  group_by(Ligand) %>%
  summarise(
    Min      = round(min(Energy, na.rm = TRUE), 3),
    Q1       = round(quantile(Energy, 0.25, na.rm = TRUE), 3),
    Median   = round(median(Energy, na.rm = TRUE), 3),
    Q3       = round(quantile(Energy, 0.75, na.rm = TRUE), 3),
    Max      = round(max(Energy, na.rm = TRUE), 3),
    Mean     = sign(mean(Energy)) * exp(mean(log(abs(Energy)))),  # geometric-mean analogue as in your code
    SD       = round(sd(Energy, na.rm = TRUE), 3),
    Variance = round(var(Energy, na.rm = TRUE), 3),
    Range    = round(max(Energy, na.rm = TRUE) - min(Energy, na.rm = TRUE), 3),
    IQR      = round(IQR(Energy, na.rm = TRUE), 3),
    Skewness = round(e1071::skewness(Energy, na.rm = TRUE), 3),
    Kurtosis = round(e1071::kurtosis(Energy, na.rm = TRUE), 3),
    .groups  = "drop"
  )

stat_controls_per_frame_SP <- controls_data_SP %>%
  group_by(Frame, Ligand) %>%
  summarise(
    Min      = round(min(Energy, na.rm = TRUE), 3),
    Q1       = round(quantile(Energy, 0.25, na.rm = TRUE), 3),
    Median   = round(median(Energy, na.rm = TRUE), 3),
    Q3       = round(quantile(Energy, 0.75, na.rm = TRUE), 3),
    Max      = round(max(Energy, na.rm = TRUE), 3),
    Mean     = sign(mean(Energy)) * exp(mean(log(abs(Energy)))),  # geometric-mean analogue
    SD       = round(sd(Energy, na.rm = TRUE), 3),
    Variance = round(var(Energy, na.rm = TRUE), 3),
    Range    = round(max(Energy, na.rm = TRUE) - min(Energy, na.rm = TRUE), 3),
    IQR      = round(IQR(Energy, na.rm = TRUE), 3),
    Skewness = round(e1071::skewness(Energy, na.rm = TRUE), 3),
    Kurtosis = round(e1071::kurtosis(Energy, na.rm = TRUE), 3),
    .groups  = "drop"
  )

# -----------------------------
# 7) Cutoffs table (DJM, AHA, BME)
# -----------------------------
cutoffs  <- c(-6.747508, -5.093705, -3.522910)
names_co <- c("DJM", "AHA", "BME")

total_molecules_SP <- geo_mean_SP %>%
  summarise(total = n_distinct(Ligand)) %>% pull(total)

tabla_resultado <- geo_mean_SP %>%
  summarise(
    `Cutoff Geo Mean` = c("without cutoff", names_co),
    `# Molecules`     = c(total_molecules_SP,
                          sum(Geometric_mean <= cutoffs[1]),
                          sum(Geometric_mean <= cutoffs[2]),
                          sum(Geometric_mean <= cutoffs[3]))
  )

# -----------------------------
# 8) Selected by AHA cutoff (optional additional sheet)
# -----------------------------
SP_selected <- geo_mean_SP %>%
  filter(Geometric_mean <= -5.093705) %>%
  select(Ligand, Geometric_mean)

# -----------------------------
# 9) Histogram plot (object renamed to 'Histogram'; aesthetics unchanged)
# -----------------------------
cutoff <- -5.093705

Histogram <- ggplot(geo_mean_SP, aes(x = Geometric_mean)) +
  geom_histogram(
    aes(fill = after_stat(x > cutoff)),
    position = "dodge",
    colour = "black",
    alpha = 1,
    bins = 50
  ) +
  scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "white")) +
  labs(x = "Docking energy \nwith SP (kcal/mol)", y = "Frequency") +
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
  coord_flip() +
  scale_y_continuous(labels = comma)

# Preview (optional)
print(Histogram)

ggsave(
  filename = file.path(script_dir, "Histogram_SP.png"),
  plot = Histogram,
  width = 7,
  height = 5,
  dpi = 600,
  units = "in",
  bg = "white"
)

# -----------------------------
# 10) Per-frame plot (object renamed to 'PerFrame'; aesthetics unchanged)
# -----------------------------
geo_mean_per_frame_SP_clean <- geo_mean_per_frame_SP %>%
  filter(is.finite(Geometric_mean)) %>%
  mutate(Frame = as.numeric(Frame)) %>%
  arrange(Frame) %>%
  mutate(Frame = factor(Frame, levels = unique(Frame)))

controls_df <- stat_controls_per_frame_SP %>%
  filter(Ligand %in% c("AHA", "BME", "DJM")) %>%
  mutate(Frame = as.numeric(Frame)) %>%
  arrange(Frame) %>%
  mutate(Frame = factor(Frame, levels = unique(geo_mean_per_frame_SP_clean$Frame))) %>%
  rename(Control = Ligand, ControlMean = Mean)

PerFrame <- ggplot(geo_mean_per_frame_SP_clean, aes(y = Frame, x = Geometric_mean)) +
  geom_boxplot(
    fill = "white",
    colour = "black",
    outlier.colour = "gray40",
    outlier.alpha = 0.3,
    outlier.size = 1.2,
    width = 0.7
  ) +
  geom_point(
    data = controls_df,
    aes(x = ControlMean, y = Frame, fill = Control),
    shape = 21, size = 2.2, colour = "black"
  ) +
  scale_fill_manual(
    name = "Controls",
    values = c(AHA = "blue", BME = "green", DJM = "red")
  ) +
  labs(x = "Docking energy with SP (kcal/mol)", y = "Frame") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 1),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title   = element_text(color = "black", face = "bold", size = 22),
    axis.text.x  = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y  = element_text(size = 16, color = "black", face = "bold"),
    axis.ticks   = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    strip.text   = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold", color = "black"),
    legend.text  = element_text(size = 16, face = "bold", color = "black"),
    legend.key.size = unit(1.2, "lines"),
    legend.spacing.x = unit(0.5, "cm")
  )

# Preview (optional)
print(PerFrame)

ggsave(
  filename = file.path(script_dir, "PerFrame_SP.png"),
  plot = PerFrame,
  width = 8,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)

# -----------------------------
# 11) Write all stats to Excel (no date suffix)
# -----------------------------
write_xlsx(
  x = list(
    stat_vs_SP                 = stat_vs_SP,
    stat_vs_per_frame_SP       = stat_vs_per_frame_SP,
    stat_controls_SP           = stat_controls_SP,
    stat_controls_per_frame_SP = stat_controls_per_frame_SP,
    tabla_resultado            = tabla_resultado,
    SP_selected                = SP_selected
  ),
  path = file.path(script_dir, "Stats_SP.xlsx")
)

# -----------------------------
# 12) Done
# -----------------------------
message("Saved plot:  ", file.path(script_dir, "Histogram.png"))
message("Saved plot:  ", file.path(script_dir, "PerFrame.png"))
message("Saved stats: ", file.path(script_dir, "Stats_SP.xlsx"))

