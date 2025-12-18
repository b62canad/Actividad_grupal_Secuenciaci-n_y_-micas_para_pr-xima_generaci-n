############################################################
# 0. Paquetes y entorno
############################################################

# Estos los instalas solo si faltan (desde CRAN)
if (!requireNamespace("readr", quietly = TRUE))
  install.packages("readr")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!requireNamespace("ggrepel", quietly = TRUE))
  install.packages("ggrepel")

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!requireNamespace("enrichR", quietly = TRUE))
  install.packages("enrichR")

# Este de Bioconductor
BiocManager::install("tximport")

# Estos vienen del entorno conda/bioconda
library(tximport)
library(DESeq2)

# Estos son de CRAN
library(readr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(enrichR)

# Limpiamos el entorno y elegimos nuestro directorio de trabajo
rm(list=ls())

setwd("C:/Users/Secuenciación/mubio03_act2/TallerGrupal_Ficheros")

############################################################
# 1. Importar cuantificaciones de Salmon
############################################################

dir_quant <- "C:/Users/Secuenciación/mubio03_act2/TallerGrupal_Ficheros/Salmon"

# Carpetas con *_quant
samples <- list.files(dir_quant, pattern = "_quant$", full.names = FALSE)

# Rutas a quant.sf
files <- file.path(dir_quant, samples, "quant.sf")

# Quitar el sufijo "_quant" para tener nombres limpios de muestra
names(files) <- sub("_quant$", "", samples)

# Comprobar que están todas
files

# Se lee el archivo del diseño experimental
sample_info = read_csv("Design.csv", locale = locale(encoding = "UTF-8"))

# Se lee el archivo de correspondencia entre transcritos y genes
tx2gene = read_tsv("Transcrito_a_Gen.tsv", col_names = FALSE)

# Asignamos nombres a las columnas del data frame
colnames(tx2gene) = c("TXNAME", "GENEID")

# Importar los resultados de cuantificación generados por Salmon
txi = tximport(files, type = "salmon", tx2gene = tx2gene)


############################################################
# 2. Definir diseño experimental (Normopeso vs Obeso2)
############################################################

# Simplificamos el nombre de la condicion Sobrepeso/Obeso2 a solamente Obeso2
sample_info$Condition[sample_info$Condition == "Sobrepeso/Obeso2"] <- "Obeso2"

# Filtrar las dos condiciones de interés
sample_info_sub <- sample_info[sample_info$Condition %in% c("Obeso2", "Normopeso"), ]

rownames(sample_info_sub) <- sample_info_sub$Sample

# Asegurar coincidencia con los nombres de columnas
common_samples <- intersect(rownames(sample_info_sub), colnames(txi$counts))

sample_info_sub <- sample_info_sub[common_samples, ]

# Redondear para DESeq2
counts_sub <- round(txi$counts[, common_samples])  

############################################################
# 3. Crear objeto DESeq2 y correr DESeq
############################################################

# Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = sample_info_sub,
  design = ~ Condition
)

# Factorizar la condición
dds$Condition <- factor(dds$Condition, levels = c("Normopeso", "Obeso2"))

# Filtrar genes de baja expresión
dds <- dds[rowSums(counts(dds)) > 10, ]

# Ejecutar DESeq2
dds <- DESeq(dds)

# Resultados Obeso2 vs Normopeso
res <- results(dds, contrast = c("Condition", "Obeso2", "Normopeso"))


# Resultados DESeq

# Muestra un resumen estadístico de los resultados
summary(res)

# Convertir el objeto res de DESeq2) a un data frame para facilitar su visualización
res_df <- as.data.frame(res)

# Ver los 20 genes más importantes ordenados por padj
res <- res[order(res$padj), ]
head(res_df, 20)

############################################################
# 4. MA plot
############################################################

# Para las anotaciones de los genes
res_df$gene <- rownames(res_df)

# Clasificación de significancia
res_df$signif <- ifelse(res_df$padj < 0.05, "Significativo", "No significativo")

ggplot(res_df, aes(x = log2(baseMean), y = log2FoldChange, color = signif)) +
  geom_point(alpha = 0.5) +
  
  # Anotacion de los genes significativos
  geom_text_repel(
    data = subset(res_df, padj < 0.05),
    aes(label = gene),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = c(
    "Significativo" = "red",
    "No significativo" = "black"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "log2(baseMean)",
    y = "log2FoldChange",
    title = "MA plot: Obeso2 vs Normopeso",
    color = "Estado"
  ) +
  theme_minimal()

############################################################
# 5. Volcano plot
############################################################

# Clasificación de significancia
res_df$signif <- "No significativo"
res_df$signif[res_df$padj < 0.05 & res_df$log2FoldChange > 1]  <- "Sobreexpresado"
res_df$signif[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Subexpresado"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  # Anotación de los genes significativos
  geom_text_repel(
    data = subset(res_df,padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = gene),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = c(
    "Subexpresado" = "royalblue",
    "No significativo" = "grey50",
    "Sobreexpresado" = "firebrick"
  )) +
  labs(
    x = "log2(Fold Change) Obeso2 vs Normopeso",
    y = "-log10(p-valor ajustado)",
    title = "Volcano plot: Obeso2 vs Normopeso",
    color = "Expresión"
  ) +
  theme_minimal()

############################################################
# 6. PCA (usando varianceStabilizingTransformation)
############################################################

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extraer datos de PCA sin dibujar (evita el error de viewports)
pcaData <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  labs(title = "PCA de muestras\n(Normopeso vs Obeso2)") +
  theme_minimal()


############################################################
# 7. Heatmap de los genes más variables
############################################################

# Creamos data frame para las anotaciones
annota <- data.frame(row.names = sample_info_sub$Sample, Condition = factor(sample_info_sub$Condition))

# Seleccionamos los 40 genes con mayor varianza
topVar <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)

mat <- assay(vsd)[topVar, ]
mat <- mat - rowMeans(mat)   # centrado por gen
pheatmap(
  mat,
  annotation_col = annota,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  scale          = "row",
  fontsize_row   = 8,
  fontsize_col   = 8,
  angle_col      = 315,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  main           = "Heatmap de genes más variables"
)


############################################################
# Enriquecimiento
############################################################

# Se filtran los genes diferencialmente expresados significativos
res_sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
genes <- rownames(res_sig)

# Se define la base de datos de enriquecimiento a utilizar, en este caso, procesos biológicos de Gene Ontology (GO)
dbs_use <- c("GO_Biological_Process_2023")

# Se realiza el análisis de enriquecimiento funcional con enrichR
enrich_res <- enrichr(genes, dbs_use)

# Extraemos los resultados obtenidos de GO Biological Process
go_bp <- enrich_res$GO_Biological_Process_2023

# Se filtran los términos GO significativamente enriquecidos
go_bp_sig <- go_bp[go_bp$Adjusted.P.value < 0.05, ]

# Visualizamos las primeras filas
head(go_bp_sig)

# Se grafican los 5 primeros resultados
top_terms <- go_bp_sig[1:5, ]
ggplot(top_terms, aes(
  x = Combined.Score,
  y = reorder(Term, Combined.Score)
)) +
  geom_point(size = 3) +
  labs(
    x = "Combined Score",
    y = "GO Biological Process"
  ) +
  theme_minimal()
