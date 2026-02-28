## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(colorRamp2)


## Load Aging Data -------------------------------------------------------------

aging_tpm <- read.csv("input_data\\aging_133_tpm.csv", 
                      row.names = 1, 
                      sep = ";") ## Import curated data from GSE113957

metadata <- read.csv("input_data\\sample_metadata_correct.csv", sep = ";")

summary(metadata$Sample == colnames(aging_tpm)) ## False, needs to be corrected

aging_tpm <- aging_tpm[,metadata$Sample] 

summary(metadata$Sample == colnames(aging_tpm)) ## Reorders counts for DESeq2

metadata <- metadata %>% 
  dplyr::mutate(age_group = case_when(
    Age %in% c(1:20) ~ "1-20y",
    Age %in% c(21:40) ~ "21-40y",
    Age %in% c(41:60) ~ "41-60y",
    Age %in% c(61:80) ~ "61-80y",
    Age %in% c(81:120) ~ "80y+",
  ))

metadata_filt <- metadata[,c(1,7)]

aging_tpm_t <- t(aging_tpm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% 
  left_join(metadata_filt, join_by(Sample))


## Get mean TPM, for each gene, for each group.
aging_tpm_t_median <- aging_tpm_t %>% 
  dplyr::group_by(age_group) %>% 
  summarise_at(2:39377, ~ round(mean(.), 2)) 

aging_tpm_t_median_scaled <- aging_tpm_t_median %>% 
  column_to_rownames(var = "age_group")

aging_tpm_t_median_scaled_true <- scale(aging_tpm_t_median_scaled)

aging_tpm_t_median_scaled_true <- as.data.frame(
  aging_tpm_t_median_scaled_true
)


## FOXM1 DNA Repair Sig. -------------------------------------------------------


foxm1_dna_repair <- read.csv("foxm1_dna_repair_net.csv", sep = ";")

foxm1_dna_repair_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = foxm1_dna_repair$gene_symbol, 
                                                 columns="ENTREZID", 
                                                 keytype="SYMBOL")

aging_tpm_heatmap <- aging_tpm_t_median_scaled_true[,foxm1_dna_repair_entrez$ENTREZID]

colnames(aging_tpm_heatmap) <- foxm1_dna_repair_entrez$SYMBOL

aging_tpm_heatmap <- aging_tpm_heatmap %>% 
  select_if(~ !any(is.na(.)))


## Plotting the Heatmap:

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2)) ## Necessary to assign values to colors

Heatmap(aging_tpm_heatmap,
        col = col_fun,
        cluster_rows = F,
        column_names_gp = gpar(fontsize = 10),
        border_gp = gpar(col = "black", lty = 1),
        heatmap_legend_param = list(
          at = c(-2, 0, 2),
          labels = c("-2", "0", "2"),
          title = "Row Z-Score",
          border = "black",
          title_position = "leftcenter-rot"))
