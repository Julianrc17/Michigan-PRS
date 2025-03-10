# Datos de Score
### Cargar el archivo en R
df <- read.table("scores.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

### Ver las primeras filas para inspeccionar
head(df)

# Datos Fenotípicos

feno <- read.table("pheno.txt", header = TRUE)
head(feno)

### Juntamos el FID y el IID ya que Plink los suele separar entonces los juntamos ya que es el nombre original de las muestras.

library(dplyr)
library(tidyr)

feno <- feno %>%
  unite("FID_IID", FID, IID, sep = "_", remove = FALSE)

feno <- feno %>%
  select(-FID, -IID)

head(feno)

### Añadimos las columnas VACCINE BATCH y HIV (FID_IID al ser la misma que sample no se añade) a nuestros datos de score o puntuación para poder realizar un análisis LIMMA posteriormente.

df_merged <- merge(df, feno[, c("FID_IID", "VACCINE", "BATCH", "HIV")], 
            by.x = "sample", by.y = "FID_IID", all.x = TRUE)
head(df_merged)
colnames(df_merged)

df_merged$HIV
df_merged$VACCINE
df_merged$BATCH

## Análisis Global Pacientes

library(dplyr)
library(broom)

prs_columns <- grep("^PGS", names(df_merged), value = TRUE)

results_list <- list()

for (prs_column in prs_columns) {
  # Crear la fórmula del modelo
  formula <- as.formula(paste(prs_column, "~ HIV + VACCINE + BATCH"))
  
  model <- lm(formula, data = df_merged)
  
  results_list[[prs_column]] <- tidy(model)
}

results_df <- bind_rows(results_list, .id = "PRS")

results_df$p.adjusted <- p.adjust(results_df$p.value, method = "BH")
head(results_df)

library(openxlsx)
output_file <- "model_results.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "Results")

writeData(wb, sheet = "Results", x = results_df, startRow = 1, startCol = 1, rowNames = FALSE)

saveWorkbook(wb, file = output_file, overwrite = TRUE)

NONVAXGEN <- read.table("nonVAXGEN_IDs.txt", header = TRUE)
head(NONVAXGEN)

## Análisis de los pacientes NONVAXGEN

df_samplename <- merge(df, feno[, c("FID_IID", "VACCINE", "BATCH", "HIV", "SampleName")], 
            by.x = "sample", by.y = "FID_IID", all.x = TRUE)

library(dplyr)

# Filtrar el conjunto de datos para excluir VAXGEN
feno_nonvaxgen <- df_samplename %>% filter(!grepl("^VAXGEN", SampleName))

head(feno_nonvaxgen)

# Obtener los nombres de todas las columnas de scores de PRS
prs_columns <- grep("^PGS", names(feno_nonvaxgen), value = TRUE)

# Lista para almacenar los resultados
results_list <- list()

# Bucle para ajustar modelos lineales para cada columna de PRS
for (prs_column in prs_columns) {
  # Crear la fórmula del modelo
  formula <- as.formula(paste(prs_column, "~ HIV + BATCH"))
  
  # Ajustar el modelo lineal
  model <- lm(formula, data = feno_nonvaxgen)
  
  # Guardar los resultados en una lista
  results_list[[prs_column]] <- tidy(model)
}

# Convertir la lista de resultados en un dataframe largo
results_nonvaxgen <- bind_rows(results_list, .id = "PRS")

# Ajustar los valores p por Benjamini-Hochberg
results_nonvaxgen$p.adjusted <- p.adjust(results_nonvaxgen$p.value, method = "BH")

# Ver los primeros resultados
head(results_nonvaxgen)



output_file <- "model_results_nonvaxgen.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "Results")
writeData(wb, sheet = "Results", x = results_nonvaxgen, startRow = 1, startCol = 1, rowNames = FALSE)
saveWorkbook(wb, file = output_file, overwrite = TRUE)

## Análisis de los pacientes VAXGEN

feno_vaxgen <- df_samplename %>% filter(grepl("^VAXGEN", SampleName))
# Obtener los nombres de todas las columnas de scores de PRS
prs_columns <- grep("^PGS", names(feno_vaxgen), value = TRUE)

# Lista para almacenar los resultados
results_list <- list()

# Bucle para ajustar modelos lineales para cada columna de PRS
for (prs_column in prs_columns) {
  # Crear la fórmula del modelo
  formula <- as.formula(paste(prs_column, "~ HIV + VACCINE + BATCH"))
  
  # Ajustar el modelo lineal
  model <- lm(formula, data = feno_vaxgen)
  
  # Guardar los resultados en una lista
  results_list[[prs_column]] <- tidy(model)
}

# Convertir la lista de resultados en un dataframe largo
results_vaxgen <- bind_rows(results_list, .id = "PRS")

results_vaxgen$p.adjusted <- p.adjust(results_vaxgen$p.value, method = "BH")

head(results_vaxgen)

output_file <- "model_results_vaxgen.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "Results")
writeData(wb, sheet = "Results", x = results_vaxgen, startRow = 1, startCol = 1, rowNames = FALSE)
saveWorkbook(wb, file = output_file, overwrite = TRUE)

## Análisis de Correlación de los Scores

library(Hmisc, lib = "/home/julian/R/library")
library(PerformanceAnalytics, lib = "/home/julian/R/library")
df_num <- df[, -1]
res <- rcorr(as.matrix(df_num))  
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# Formatear la matriz de correlación en tabla
corr_table <- flattenCorrMatrix(res$r, res$P)

# Ver las primeras filas
head(corr_table)
library(dplyr)

# Lista de columnas que quieres conservar
columnas_deseadas <- c("PGS003131", "PGS002046", "PGS003130", "PGS001887", "PGS002802", "PGS001788", 
                       "PGS003835", "PGS002100", "PGS000812", "PGS002801", "PGS000659", "PGS002805", 
                       "PGS003490", "PGS003641", "PGS003519", "PGS003152", "PGS003888", "PGS002957", 
                       "PGS004065", "PGS002804", "PGS001241", "PGS001688", "PGS003147", "PGS002858", 
                       "PGS003791", "PGS000026", "PGS002962", "PGS003836", "PGS000891", "PGS001783", 
                       "PGS002333", "PGS001168", "PGS000876", "PGS000312", "PGS000989", "PGS001513", 
                       "PGS000479", "PGS002900", "PGS000718", "PGS001160", "PGS002904", "PGS001763", 
                       "PGS002524", "PGS000193", "PGS002903", "PGS000988", "PGS001351", "PGS001338", 
                       "PGS002699", "PGS001155", "PGS000334", "PGS001063", "PGS003125", "PGS001002", 
                       "PGS003479", "PGS001392", "PGS002899", "PGS002457", "PGS002672")

# Filtrar solo las columnas deseadas
df_filtrado <- df_merged %>% select(any_of(columnas_deseadas))

# Ver las primeras filas del nuevo dataframe
head(df_filtrado)
library(corrplot, lib = "/home/julian/R/library")

cor_matrix <- cor(df_filtrado, use = "pairwise.complete.obs", method = "pearson")

# Visualizar la matriz de correlación con un heatmap
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.6, tl.col = "black", col = colorRampPalette(c("blue", "white", "red"))(200))
