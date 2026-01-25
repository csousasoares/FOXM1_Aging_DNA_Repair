library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(forcats)
library(limma)
library(msigdbr)
library(org.Hs.eg.db)
library(stringr)
library(patchwork)
library(tidytree)
library(ggrastr)

set.seed(123)

mm_BP_sets <- msigdbr(
  species = "Homo sapiens")
mm_BP_sets

msigdbr_t2g = mm_BP_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_t2g

misgdbr_t2g_filtered <- filter(msigdbr_t2g, gs_name %in% c("GOBP_BASE_EXCISION_REPAIR",
                                                           "GOBP_DNA_REPAIR",
                                                           "GOBP_DNA_SYNTHESIS_INVOLVED_IN_DNA_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
                                                           "GOBP_INTERSTRAND_CROSS_LINK_REPAIR",
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
                                                           "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "REACTOME_DNA_REPAIR",
                                                           "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
                                                           "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
                                                           "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS",
                                                           "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
                                                           "WP_BASE_EXCISION_REPAIR",
                                                           "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                                                           "WP_NUCLEOTIDE_EXCISION_REPAIR"))


aging_133 <- read.table("aging_ranked_list.txt", sep="\t", header=T, row.names = 1)
head(aging_133)  
gene_list_aging_133 <- aging_133$log2fc
names(gene_list_aging_133) <- rownames(aging_133)
gene_list_aging_133


go_analysis_aging_133 <- GSEA(gene_list_aging_133,
                               exponent = 1,
                               minGSSize = 0,
                               maxGSSize = 1000,
                               eps = 1e-300,
                               pvalueCutoff = 1,
                               pAdjustMethod = "fdr",
                               gson = NULL,
                               TERM2GENE = misgdbr_t2g_filtered,
                               verbose = T,
                               seed = F,
                               by = "fgsea")

df <- as.data.frame(go_analysis_aging_133)
dotplot(go_analysis_aging_133, showCategory=40, 
        title = "Over 40y vs Under 40y DNA Repair", 
        orderBy = "NES", x = "NES", color = "p.adjust", 
        font.size = 7, label_format = 60) + facet_grid(.~.sign)
write.csv(go_analysis_aging_133, "output_data\\\go_analysis_aging_133_DNA_REPAIR.csv", 
          row.names = F)
gseaplot2(go_analysis_aging_133, geneSetID = 4,
          title = go_analysis_aging_133$Description[4],
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.5,0.5),
          subplots = 1:3,
          pvalue_table = T,
          ES_geom = "line")

##HALMARKS

mm_BP_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H")
mm_BP_sets

msigdbr_t2g = mm_BP_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_t2g

go_analysis_aging_133_hall <- GSEA(gene_list_aging_133,
                              exponent = 1,
                              minGSSize = 0,
                              maxGSSize = 1000,
                              eps = 1e-300,
                              pvalueCutoff = 1,
                              pAdjustMethod = "fdr",
                              gson = NULL,
                              TERM2GENE = msigdbr_t2g,
                              verbose = T,
                              seed = F,
                              by = "fgsea",
                              nPermSimple = 10000)
df2 <- as.data.frame(go_analysis_aging_133_hall)
dotplot(go_analysis_aging_133_hall, showCategory=40, title = "Over 40y vs Under 40y HALL", orderBy = "NES", x = "NES", color = "p.adjust", font.size = 7, label_format = 60) + facet_grid(.~.sign)
write.csv(go_analysis_aging_133_hall, "output_data\\go_analysis_aging_133_HALL.csv", row.names = F)
gseaplot2(go_analysis_aging_133_hall, geneSetID = 1,
          title = go_analysis_aging_133_hall$Description[1],
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.5,0.5),
          subplots = 1:3,
          pvalue_table = T,
          ES_geom = "line")

##Load GSEA RData

trace("gseaplot2", edit = T)

dna_rep_plot <- gseaplot2(go_analysis_aging_133, geneSetID = c(1,4,7,9,16),
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.3,0.5),
          subplots = 1:2,
          pvalue_table = F,
          ES_geom = "line")

dna_rep_plot
str(dna_rep_plot$layers)


ggrastr::rasterise(dna_rep_plot, layers = "Line", dpi = 100)



dsb_rep_plot <- gseaplot2(go_analysis_aging_133, geneSetID = c(3,5),
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.3,0.5),
          subplots = 1:2,
          pvalue_table = F,
          ES_geom = "line")

hr_nhej_plot <- gseaplot2(go_analysis_aging_133, geneSetID = c(2, 19, 25),
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.3,0.5),
          subplots = 1:2,
          pvalue_table = T,
          ES_geom = "line")

hr_nhej_plot

dsb_hr_nhej_plot <- gseaplot2(go_analysis_aging_133, geneSetID = c(2, 3,5,19, 25),
                          base_size = 12,
                          color = "green",
                          rel_heights = c(2,0.3,0.5),
                          subplots = 1:2,
                          pvalue_table = F,
                          ES_geom = "line")

dsb_hr_nhej_plot


dsb_hr_nhej_plot

ber_plot <- gseaplot2(go_analysis_aging_133, geneSetID = c(8, 10, 11, 13),
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.3,0.5),
          subplots = 1:2,
          pvalue_table = F,
          ES_geom = "line")

ber_plot



ner_plot <- gseaplot2(go_analysis_aging_133, geneSetID = c(15,23,24, 30),
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.3,0.5),
          subplots = 1:2,
          pvalue_table = F,
          ES_geom = "line")

ner_plot


icl_plot <- gseaplot2(go_analysis_aging_133, geneSetID = c(17),
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.3,0.5),
          subplots = 1:2,
          pvalue_table = F,
          ES_geom = "line")


icl_plot

library(patchwork)

# Combine the first list into a single plot
ner_plot_combined <- wrap_plots(ner_plot, ncol = 1)

# Combine the second list into a single plot
icl_plot_combined <- wrap_plots(icl_plot, ncol = 1)

# Now combine both vertically or horizontally
final_plot <- ner_plot_combined | icl_plot_combined  # side by side
# Or use / for stacking vertically: ner_plot_combined / icl_plot_combined

# Print
final_plot
