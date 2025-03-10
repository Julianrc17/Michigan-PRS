# Score Data Analysis in R

## Overview
This project analyzes score data using R, integrating phenotypic data and performing statistical analyses, including linear modeling and correlation analysis.

## Loading Score Data
```r
# Load the score data
df <- read.table("scores.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
# View first rows
head(df)
```

## Loading Phenotypic Data
```r
feno <- read.table("pheno.txt", header = TRUE)
head(feno)
```

## Merging FID and IID
PLINK separates FID and IID; thus, we merge them.
```r
library(dplyr)
library(tidyr)

feno <- feno %>%
  unite("FID_IID", FID, IID, sep = "_", remove = FALSE) %>%
  select(-FID, -IID)

head(feno)
```

## Merging Phenotypic Data with Score Data
```r
df_merged <- merge(df, feno[, c("FID_IID", "VACCINE", "BATCH", "HIV")],
            by.x = "sample", by.y = "FID_IID", all.x = TRUE)
head(df_merged)
```

## Global Patient Analysis
```r
library(broom)
prs_columns <- grep("^PGS", names(df_merged), value = TRUE)
results_list <- list()

for (prs_column in prs_columns) {
  formula <- as.formula(paste(prs_column, "~ HIV + VACCINE + BATCH"))
  model <- lm(formula, data = df_merged)
  results_list[[prs_column]] <- tidy(model)
}

results_df <- bind_rows(results_list, .id = "PRS")
results_df$p.adjusted <- p.adjust(results_df$p.value, method = "BH")
```

### Save Results
```r
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Results")
writeData(wb, sheet = "Results", x = results_df, startRow = 1, startCol = 1)
saveWorkbook(wb, file = "model_results.xlsx", overwrite = TRUE)
```

## NONVAXGEN Patient Analysis
```r
NONVAXGEN <- read.table("nonVAXGEN_IDs.txt", header = TRUE)
feno_nonvaxgen <- df_samplename %>% filter(!grepl("^VAXGEN", SampleName))
```

## VAXGEN Patient Analysis
```r
feno_vaxgen <- df_samplename %>% filter(grepl("^VAXGEN", SampleName))
```

## Score Correlation Analysis
```r
library(Hmisc)
res <- rcorr(as.matrix(df[, -1]))
```

### Format Correlation Matrix
```r
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = cormat[ut],
    p = pmat[ut]
  )
}
corr_table <- flattenCorrMatrix(res$r, res$P)
```

### Filter Selected PRS Columns
```r
columnas_deseadas <- c("PGS003131", "PGS002046", "PGS003130", "PGS001887", "PGS002802")
df_filtrado <- df_merged %>% select(any_of(columnas_deseadas))
head(df_filtrado)
```

