# -------------------------------
# Load required libraries
# -------------------------------
library(dplyr)
library(ggplot2)
library(scales)


# -------------------------------
# Load the dataset from an .RDS file
# -------------------------------
# Make sure the .rds file is in your working directory or provide the full path
datos_long <- readRDS("MegaDataML.rds")


# -------------------------------
# Create a bar plot comparing model predictions
# -------------------------------

# Summarize counts and percentages of predictions per model
summary_models <- datos_long %>%
  group_by(Model, Prediction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Model) %>%
  mutate(Percentage = Count / sum(Count) * 100)

print(summary_models)

# -------------------------------
# Generate a comparative bar plot
# -------------------------------
comparison_models <- ggplot(summary_models, aes(x = Model, y = Percentage, fill = Prediction)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 4,
    family = "Arial",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("yes" = "#2E8B57", "no" = "#DC143C"),
    labels = c("yes" = "UI", "no" = "nUI")
  ) +
  labs(
    y = "Percentage (%)",
    x = "ML Model",
    fill = "Prediction"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    axis.text.x  = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y  = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text  = element_text(size = 11),
    axis.ticks   = element_line(colour = "black", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.position = "top"
  )

# Display the plot
comparison_models

# Save the plot 
ggsave(
  filename = "comparison_models.png",
  plot = comparison_models,
  width = 7,
  height = 5,
  dpi = 600,
  units = "in",
  bg = "white"
)
