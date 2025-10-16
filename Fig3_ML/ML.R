library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


files <- list.files(path = "~/Documentos/MLVS_Triazols/CodigosML/Excels", pattern = "*.xlsx", full.names = T)

data <- lapply(files, read_excel) %>% bind_rows()

data_long <- data %>%
  pivot_longer(cols = starts_with("dt_"):starts_with("rf_"),
               names_to = "modelo", values_to = "prediccion")

####################### contar los yes y los no por modelo
conteo <- data_long %>%
  count(modelo, prediccion)

ggplot(conteo, aes(x = modelo, y = n, fill = prediccion)) +
  geom_col(position = "stack") +
  labs(y = "Número de ligandos", x = "Modelo", fill = "Predicción") +
  theme_minimal()



####################################################################################################################################################
#barplot

columnas_modelos <- c("dt_boruta_5_50", "xgb_nfs_25_50", "dt_xgb_5", "knn_xgb_5_50", "rf_boruta_10_50")

datos_largo <- data%>%
  pivot_longer(cols = all_of(columnas_modelos),
               names_to = "Modelo", 
               values_to = "Predicción") %>%
  mutate(Predicción = factor(Predicción, levels = c("yes", "no")))

resumen_modelos <- datos_largo %>%
  group_by(Modelo, Predicción) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Modelo) %>%
  mutate(Percentage = Count/sum(Count)*100)

print(resumen_modelos)

#opción 1
library(ggplot2)
library(scales)

comparacion_modelos <- ggplot(resumen_modelos, aes(x = Modelo, y = Percentage, fill = Predicción)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 4, family = "Arial", fontface = "bold") +
  scale_fill_manual(values = c("yes" = "#2E8B57", "no" = "#DC143C"),
                    labels = c("yes" = "UI", "no" = "nUI")) +
  labs(
    #title = "Comparación de Predicciones entre Modelos",
    #subtitle = "Distribución de predicciones Yes/No por modelo",
    y = "Percentage (%)", 
    x = "ML model",
    fill = "Prediction"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_classic(base_family = "Arial") +
  theme(
    #plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    #plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    axis.text.x  = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y  = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
    legend.position = "top"
  )

comparacion_modelos
ggsave(
  filename = "comparacion_modelos.png",  # nombre del archivo
  plot = comparacion_modelos,            # objeto del gráfico
  width = 7,                             # ancho (igual que hypo_db)
  height = 5,                            # alto (igual que hypo_db)
  dpi = 600,                             # resolución para journals (ACS)
  units = "in",                          # pulgadas
  bg = "white"                           # fondo blanco
)


#estadísticas resumen

cat("=== REPORTE DE ANÁLISIS DE PREDICCIONES ===\n")
cat("Total de ligandos analizados:", nrow(data), "\n")
cat("Número de archivos:", length(files), "\n")
cat("Modelos incluidos:", paste(columnas_modelos, collapse = ", "), "\n\n")

# Estadísticas de consenso
unanime_yes <- sum(datos_consenso$votos_yes == 5)
unanime_no <- sum(datos_consenso$votos_no == 5)
mayoria_yes <- sum(datos_consenso$consenso == "Mayoría Yes")
mayoria_no <- sum(datos_consenso$consenso == "Mayoría No")

cat("--- ESTADÍSTICAS DE CONSENSO ---\n")
cat("Acuerdo unánime 'Yes':", unanime_yes, 
    paste0("(", round(unanime_yes/nrow(data)*100, 1), "%)\n"))
cat("Acuerdo unánime 'No':", unanime_no, 
    paste0("(", round(unanime_no/nrow(data)*100, 1), "%)\n"))
cat("Mayoría 'Yes':", mayoria_yes, 
    paste0("(", round(mayoria_yes/nrow(data)*100, 1), "%)\n"))
cat("Mayoría 'No':", mayoria_no, 
    paste0("(", round(mayoria_no/nrow(data)*100, 1), "%)\n"))
cat("Total con consenso claro:", unanime_yes + unanime_no + mayoria_yes + mayoria_no,
    paste0("(", round((unanime_yes + unanime_no + mayoria_yes + mayoria_no)/nrow(data)*100, 1), "%)\n\n"))

# Porcentaje de yes por modelo
cat("--- TENDENCIA POR MODELO ---\n")
for(modelo in columnas_modelos) {
  perc_yes <- mean(data[[modelo]] == "yes", na.rm = TRUE) * 100
  cat(modelo, ":", round(perc_yes, 1), "% Yes\n")
}





################################ UPSET

# Preparar los datos
datos_upset <- data %>%
  mutate(across(all_of(columnas_modelos), ~ .x == "yes")) %>%
  select(all_of(columnas_modelos))


# Calcular set sizes (totales por modelo)
set_sizes_df <- datos_upset %>%
  summarise(across(all_of(columnas_modelos), sum)) %>%
  pivot_longer(everything(), names_to = "modelo", values_to = "count") %>%
  mutate(modelo = factor(modelo, levels = columnas_modelos))

plot_setsizes <- ggplot(set_sizes_df, aes(x = count, y = modelo)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7, width = 0.7) +
  geom_text(aes(label = scales::comma(count)), hjust = -0.2, size = 6, fontface = "bold") +  
  labs(
    x = "Set Size",
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 1),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title   = element_text(color = "black", face = "bold", size = 22),
    axis.text.x  = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y  = element_text(size = 16, color = "black", face = "bold"),
    legend.position = "none",
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    strip.text = element_blank(),
    plot.margin = margin(5, 10, 5, 10)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15)),
    labels = scales::comma  # Esto formatea el eje X con separadores de miles
  )

plot_setsizes

#############
# Calcular TODAS las intersecciones
intersecciones_df <- calcular_intersecciones(datos_upset, columnas_modelos)

plot_intersecciones <- ggplot(intersecciones_df, 
                              aes(x = reorder(combinacion, frecuencia), y = frecuencia)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = scales::comma(frecuencia)), vjust = -0.5, size = 3.5) +
  labs(
    y = "Frequency",
    x = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 1),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title   = element_text(color = "black", face = "bold", size = 22),
    axis.text.x  = element_blank(),  # Ocultar texto del eje X
    axis.text.y  = element_text(size = 16, color = "black", face = "bold"),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1), 
    strip.text = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    plot.margin = margin(10, 10, 15, 10)  
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.00, 0.15)),  
    labels = scales::comma
  ) +
  coord_cartesian(ylim = c(0, NA))  

plot_intersecciones

############
# Crear matriz de presencia 
crear_matriz_presencia_ordenada <- function(intersecciones_df, modelos) {
  matriz_list <- list()
  
  for(i in 1:nrow(intersecciones_df)) {
    comb_str <- as.character(intersecciones_df$combinacion[i])
    comb_modelos <- strsplit(comb_str, " & ")[[1]]
    
    for(modelo in modelos) {
      presente <- ifelse(modelo %in% comb_modelos, 1, 0)
      matriz_list[[length(matriz_list) + 1]] <- data.frame(
        combinacion = comb_str,
        modelo = modelo,
        presente = presente
      )
    }
  }
  
  matriz_df <- bind_rows(matriz_list) %>%
    mutate(
      combinacion = factor(combinacion, levels = intersecciones_df$combinacion),
      modelo = factor(modelo, levels = modelos)
    )
  
  return(matriz_df)
}

matriz_df <- crear_matriz_presencia_ordenada(intersecciones_df, columnas_modelos)

plot_matriz <- ggplot(matriz_df, aes(x = combinacion, y = modelo, color = factor(presente))) +
  geom_point(size = 4, shape = 19, alpha = 0.8) +  # Círculos sólidos
  scale_color_manual(values = c("0" = "lightgray", "1" = "black"), guide = "none") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 1),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title   = element_text(color = "black", face = "bold", size = 22),
    axis.text.x  = element_blank(),  # Ocultar texto X
    axis.text.y  = element_text(size = 16, color = "black", face = "bold"),  
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.5),  
    legend.position = "none",
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),  
    strip.text = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    plot.margin = margin(5, 10, 5, 10) 
  )

plot_matriz
##### unir todo
library(patchwork)

# Definir paleta de colores pastel para 5 modelos
paleta_pastel <- c(
  "#D6B3C9",  
  "#838091",   
  "#9BC69F",  
  "#A8C33B",  
  "#CDB586"   
)

# Asegurar que todos los gráficos tengan el mismo orden
intersecciones_df <- intersecciones_df %>%
  mutate(combinacion = factor(combinacion, levels = combinacion))

matriz_df <- matriz_df %>%
  mutate(combinacion = factor(combinacion, levels = levels(intersecciones_df$combinacion)))

set_sizes_df <- set_sizes_df %>%
  mutate(modelo = factor(modelo, levels = columnas_modelos)) %>%
  # Añadir colores según la posición del modelo
  mutate(color_modelo = paleta_pastel[as.numeric(modelo)])

# Añadir colores a la matriz_df también
matriz_df <- matriz_df %>%
  mutate(color_modelo = paleta_pastel[as.numeric(modelo)])

# 1. Gráfico de Intersecciones (Barras verticales - ARRIBA DERECHA)
plot_intersecciones <- ggplot(intersecciones_df, 
                              aes(x = combinacion, y = frecuencia)) +
  geom_bar(stat = "identity", fill = "gray", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = scales::comma(frecuencia)), vjust = -0.5, size = 3) +
  labs(y = "Frequency", x = NULL) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.margin = margin(5, 10, 0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma)

# 2. Gráfico de Set Sizes CON COLORES POR MODELO
plot_setsizes <- ggplot(set_sizes_df, aes(x = count, y = modelo, fill = modelo)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = scales::comma(count)), hjust = -0.4, size = 3.5, fontface = "bold") + 
  labs(x = "Set Size", y = NULL) +
  scale_fill_manual(values = paleta_pastel, guide = "none") +  
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black", face = "bold"),
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.margin = margin(0, 0, 5, 10)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0)),  
    labels = scales::comma,
    trans = "reverse"  
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) 

# 3. Matriz de Presencia CON COLORES POR MODELO
plot_matriz <- ggplot(matriz_df, aes(x = combinacion, y = modelo)) +
  # Líneas verticales que conectan puntos "yes"
  geom_line(
    data = matriz_df %>% filter(presente == 1),
    aes(group = combinacion),
    color = "black",  
    alpha = 0.6,
    linewidth = 0.5,
    linetype = "dashed"  
  ) +
  # Puntos con colores por modelo y tamaños diferentes según presencia
  geom_point(
    aes(color = modelo, size = factor(presente), fill = modelo),  
    alpha = 0.9,
    shape = 21  
  ) +
  scale_color_manual(values = paleta_pastel, guide = "none") +  
  scale_fill_manual(values = paleta_pastel, guide = "none") +   
  scale_size_manual(
    values = c("0" = 2.5, "1" = 4.5),  
    guide = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.ticks.y = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.margin = margin(0, 10, 5, 0)
  )

# 4. Crear un plot vacío para el espacio superior izquierdo
plot_vacio <- ggplot() + 
  theme_void() +
  theme(panel.border = element_rect(colour = "white", fill = NA, size = 0.0))

# 5. Unir todo con el layout correcto
design <- "
ABBB
CDDD
"

plot_final <- plot_vacio + plot_intersecciones + plot_setsizes + plot_matriz +
  plot_layout(design = design, 
              widths = c(1, 3),  # Izquierda 1/4, Derecha 3/4
              heights = c(1, 1.2)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

print(plot_final)

ggsave("upset_ML.png", plot_final, width = 18, height = 8, dpi = 600)
