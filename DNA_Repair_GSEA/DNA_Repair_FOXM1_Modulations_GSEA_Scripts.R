##Load necessary packages

library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

mm_BP_sets <- msigdbr(
  species = "Homo sapiens")
mm_BP_sets

head(mm_BP_sets)
colnames(mm_BP_sets)

msigdbr_t2g = mm_BP_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_t2g ##This should have loaded all MSigDB pathways from the msigdbr package as a dataframe, using gene symbols

msigdbr_t2g_filtered <- filter(msigdbr_t2g, gs_name %in% c("GOBP_BASE_EXCISION_REPAIR",
                                                           "GOBP_DNA_REPAIR",
                                                           "GOBP_DNA_SYNTHESIS_INVOLVED_IN_DNA_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
                                                           "GOBP_INTERSTRAND_CROSS_LINK_REPAIR",
                                                           "GOBP_MISMATCH_REPAIR",
                                                           "GOBP_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
                                                           "GOBP_POSTREPLICATION_REPAIR",
                                                           "GOBP_RECOMBINATIONAL_REPAIR",
                                                           "GOBP_REGULATION_OF_DNA_REPAIR",
                                                           "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
                                                           "HALLMARK_DNA_REPAIR",
                                                           "KAUFFMANN_DNA_REPAIR_GENES",
                                                           "KEGG_BASE_EXCISION_REPAIR",
                                                           "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "REACTOME_BASE_EXCISION_REPAIR",
                                                           "REACTOME_BASE_EXCISION_REPAIR_AP_SITE_FORMATION",
                                                           "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "REACTOME_DNA_REPAIR",
                                                           "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
                                                           "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
                                                           "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS",
                                                           "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
                                                           "WP_BASE_EXCISION_REPAIR",
                                                           "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                                                           "WP_NUCLEOTIDE_EXCISION_REPAIR")) ## Filter relevant DNA Repair Pathways
head(msigdbr_t2g_filtered)

##Perform GSEA for FOXM1 RNAi - also known as siFOXM1

foxm1_RNAi_new <- read.table("siFOXM1_2024_ranked.txt", sep="\t", header=T, row.names = 1)
head(foxm1_RNAi_new)  
gene_list_FOXM1_RNAi <- foxm1_RNAi_new$log2fc
names(gene_list_FOXM1_RNAi ) <- rownames(foxm1_RNAi_new)
gene_list_FOXM1_RNAi


go_analysis_FOXM1_RNAi <- GSEA(gene_list_FOXM1_RNAi,
                             exponent = 1,
                             minGSSize = 30,
                             maxGSSize = 800,
                             eps = 1e-300,
                             pvalueCutoff = 1,
                             pAdjustMethod = "fdr",
                             gson = NULL,
                             TERM2GENE = msigdbr_t2g_filtered,
                             verbose = T,
                             seed = F,
                             by = "fgsea")

df_FOXM1_RNAi <- as.data.frame(go_analysis_FOXM1_RNAi)
write.csv(go_analysis_FOXM1_RNAi, "siFOXM1_2024_DNA_Repair.csv", row.names = F)
gseaplot2(go_analysis_FOXM1_RNAi, geneSetID = 1,
          title = go_analysis_FOXM1_RNAi$Description[1],
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.5,0.5),
          subplots = 1:3,
          pvalue_table = T,
          ES_geom = "line")

dotplot(go_analysis_FOXM1_RNAi, showCategory=40, title = "FOXM1 RNAi BPs", 
        orderBy = "NES", 
        x = "NES", 
        color = "p.adjust", 
        font.size = 7, 
        label_format = 60) + facet_grid(.~.sign)


devtools::session_info()
