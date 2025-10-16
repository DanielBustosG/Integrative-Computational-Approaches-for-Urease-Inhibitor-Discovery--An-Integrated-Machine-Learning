library(readxl)
library(ggplot2)
library(dplyr)
library(flextable)
library(knitr)

#-----------------------------------------------------------------------
# 1. CONFIGURACIÓN DE ENTRADA Y SALIDA
#-----------------------------------------------------------------------

# Usa una ruta relativa para el archivo (debe estar en el mismo directorio del script)
input_file <- "rmsd_prot.xlsx"

# Carpeta de salida (creada dentro del directorio actual)
output_folder <- file.path(getwd(), "Results_RMSD", "RMSD")

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

#-----------------------------------------------------------------------
# 2. PROCESAMIENTO DE DATOS
#-----------------------------------------------------------------------
process_sheet <- function(sheet_name) {
  data <- read_excel(input_file, sheet = sheet_name)
  colnames(data) <- c("RMSD")
  df <- data.frame(Frame = 1:nrow(data), mean = data$RMSD)
  df$Protein <- sheet_name
  return(df)
}

# Listado de hojas
sheets <- c("AHA", "BME", "DJM", "CA1", "CA2", "CA3", "CA4", "CA5", "CA6", "CA7")

# Combinar todas las hojas
data_final <- do.call(rbind, lapply(sheets, process_sheet))
data_final$Protein <- factor(data_final$Protein, levels = sheets)

#-----------------------------------------------------------------------
# 4. GRÁFICO DE LÍNEAS (LINE PLOT)
#-----------------------------------------------------------------------
linePlot <- ggplot(data = data_final, aes(x = Frame, y = mean, group = Protein, color = Protein)) +
  geom_line(linewidth = 1) +
  scale_color_grey(start = 0.1, end = 0.8) +
  theme_test() +
  scale_x_continuous(name = "Frame", limits = c(0, NA)) +
  scale_y_continuous(name = "Ligand RMSD (Å)", limits = c(0, 8), breaks = seq(from = 0, to = 8, by = 2)) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", face = "bold", size = 19),
    axis.text.x = element_text(colour = "black", size = 15),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )
linePlot

# --- CORRECCIÓN AQUÍ ---
# Guardar el GRÁFICO DE LÍNEAS (objeto: linePlot, nombre: RMSD_ligand.png)
ggsave(
  filename = file.path(output_folder, "RMSD_ligand.png"),
  plot = linePlot,
  width = 10, height = 6, dpi = 300
)

#-----------------------------------------------------------------------
# 5. PREPARACIÓN DE DATOS PARA BOXPLOT Y TABLA
#-----------------------------------------------------------------------
library(readxl)

# Archivo de entrada (en el mismo directorio del script)
input_file <- "rmsd_prot.xlsx"

# Listado de hojas
sheets <- c("AHA", "BME", "DJM", "CA1", "CA2", "CA3", "CA4", "CA5", "CA6", "CA7")

# Lectura de cada hoja y combinación
sistemas <- lapply(sheets, function(s) {
  df <- read_excel(input_file, sheet = s)
  colnames(df) <- c("RMSD")
  df$sistema <- s
  return(df)
})

datos_comb <- do.call(rbind, sistemas)
datos_comb$sistema <- factor(datos_comb$sistema, levels = sheets)

#-----------------------------------------------------------------------
# 6. GRÁFICO BOXPLOT (ESTILO BLANCO Y NEGRO)
#-----------------------------------------------------------------------
boxPlot <- ggplot(datos_comb, aes(x = sistema, y = RMSD)) +
  geom_boxplot(fill = "white", color = "black", outlier.size = 1, outlier.shape = 19) +
  theme_test() +
  scale_x_discrete(name = "Candidates") +
  scale_y_continuous(name = "Ligand RMSD (Å)") +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "black", face = "bold", size = 19),
    axis.text.x = element_text(colour = "black", size = 15, angle = 0, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

# Muestra el gráfico en RStudio
boxPlot

# --- CORRECCIÓN AQUÍ ---
# Guardar el BOXPLOT (objeto: boxPlot, nombre: RMSD_ligand-boxplot.png)
ggsave(
  filename = file.path(output_folder, "RMSD_ligand-boxplot.png"),
  plot = boxPlot,
  width = 7, height = 5, dpi = 600, units = "in", bg = "white"
)


#-----------------------------------------------------------------------
# 7. TABLA DE ESTADÍSTICAS
#-----------------------------------------------------------------------
resultados <- datos_comb %>%
  group_by(sistema) %>%
  summarise(
    Q1 = quantile(RMSD, 0.25),
    Median = median(RMSD),
    Q3 = quantile(RMSD, 0.75),
    Mean = mean(RMSD),
    S.D. = sd(RMSD)
  ) %>%
  rename(Candidates = sistema)

resultados

# Crea y muestra la tabla en el Viewer
flex_table1 <- flextable(resultados) %>%
  colformat_double(j = c("Q1", "Median", "Q3", "Mean", "S.D."), digits = 3) %>%
  set_table_properties(layout = "autofit")

print(flex_table1)

# Crea y muestra la tabla en la consola
print(kable(resultados, digits = 3))