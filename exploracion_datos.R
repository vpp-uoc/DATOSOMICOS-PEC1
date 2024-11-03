library(readxl)
library(dplyr)
library(SummarizedExperiment)
file_path <- "C:/Users/Vane/Desktop/Archivos R/Datos Ómicos/TIO2+PTYR-human-MSS+MSIvsPD.XLSX"
data_matrix <- read_excel(file_path, sheet = "originalData")
metadata <- read_excel(file_path, sheet = "targets")
metadata <- as.data.frame(metadata)

# Procesamiento y Limpieza de los Datos
quantitative_data <- as.matrix(data_matrix[, 7:ncol(data_matrix)])
quantitative_data <- apply(quantitative_data, 2, function(x) as.numeric(as.character(x)))
if (any(is.na(quantitative_data))) {
  print("Existen valores NA en los datos. Se reemplazarán con 0.")
  quantitative_data[is.na(quantitative_data)] <- 0
}

rownames(quantitative_data) <- data_matrix$Accession
colnames(quantitative_data) <- metadata$Sample

# Creación del Contenedor SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = quantitative_data),
  colData = metadata,
  rowData = data_matrix[, 1:6]  # Información de las filas
)

# Exploración de los Datos
print(summary(se))
summary(assay(se))
hist(assay(se)[, 1], main = "Histograma de la primera muestra", xlab = "Intensidad", col = "blue")

library(ggplot2)
library(ggfortify)
autoplot(prcomp(t(assay(se))), data = as.data.frame(colData(se)), colour = 'Phenotype')

# Datos del grupo MSS
mss_data <- assay(se)[, colData(se)$Phenotype == "MSS"]
pd_data <- assay(se)[, colData(se)$Phenotype == "PD"]
par(mfrow = c(1, 2))  # Dividir la ventana gráfica en 2 columnas
hist(mss_data, main = "Histograma del Grupo MSS", xlab = "Intensidad", col = "blue")
hist(pd_data, main = "Histograma del Grupo PD", xlab = "Intensidad", col = "red")
boxplot(assay(se), main = "Boxplot de las Muestras", xlab = "Muestras", ylab = "Intensidad", col = "lightblue", las = 2)

# Realizamos un ANOVA para cada fosfopéptido
anova_results <- apply(assay(se), 1, function(x) {
  fit <- aov(x ~ colData(se)$Phenotype)
  summary(fit)[[1]][["Pr(>F)"]][1]
})
head(anova_results)

# Seleccionamos los fosfopéptidos con p-valor significativo y válidos
signif_peptides <- rownames(assay(se))[!is.na(anova_results) & anova_results < 0.05]
signif_peptides
if (length(signif_peptides) > 0) {
  # Subset de fosfopéptidos significativos
  signif_data <- assay(se)[signif_peptides, ]
  # heatmap
  library(pheatmap)
  pheatmap(signif_data, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
           main = "Heatmap de Fosfopéptidos Significativos")
} else {
  print("No se encontraron fosfopéptidos significativos con un p-valor menor a 0.05.")
}

# Seleccionamos los 50 fosfopéptidos más significativos para que podamos visualizar mejor el heatmap
top_signif_peptides <- head(signif_peptides, 50)
top_signif_data <- assay(se)[top_signif_peptides, ]
pheatmap(top_signif_data, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         main = "Heatmap de los 50 Fosfopéptidos Más Significativos",
         fontsize_row = 8, fontsize_col = 8, cellwidth = 10, cellheight = 10)