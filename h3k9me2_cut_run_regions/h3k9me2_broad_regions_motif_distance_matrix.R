set.seed(123)

library(tidyverse)
library(ComplexHeatmap)

test <- read.csv("test.csv", header = T, sep = ";")
motif_table <- test %>%
  distinct(TF, CONSENSUS) %>%
  filter(!is.na(CONSENSUS))

names <- test$TF

# Build a pairwise distance matrix based on string distance:
motif_dist <- stringdist::stringdistmatrix(
  motif_table$CONSENSUS, 
  method = "lv")  # Levenshtein distance

motif_dist

motif_dist <- stringdist::stringdistmatrix(
  motif_table$CONSENSUS,
  motif_table$CONSENSUS,
  method = "lv"
)

rownames(motif_dist) <- motif_table$TF
colnames(motif_dist) <- motif_table$TF


col_fun = circlize::colorRamp2(c(1, 7), c("red" ,"white"))

pval <- test$pval_dndk
pval <- ifelse(pval == 0, 0.0001, pval)
pval2 <- -log10(pval)
pval2

pval_3 <- test$pval_sifoxm1
pval_3 <- ifelse(pval_3 == 0, 0.0001, pval_3)
pval4 <- -log10(pval_3)
pval4

row_ha_2 = rowAnnotation(bar1=pval4, bar2 = pval2)

Heatmap(motif_dist, use_raster = T, raster_quality = 2,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        col = col_fun,
        right_annotation = row_ha_2)

sessionInfo()
