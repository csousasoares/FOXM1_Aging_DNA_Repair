##Load Packages ----------------------------------------------------------------

library(rGREAT)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(simplifyEnrichment)
library(magick)
library(GenomicRanges)
library(rtracklayer)
library(msigdbr)
library(dplyr)


set.seed(123)

## Load Data for rGREAT Analysis -----------------------------------------------

## Load H3K9me2 Regions Overlapping with Enhancers obtained
## Using bedtools intersect and the ENCODE cCREs Database

sifoxm1_enhancers_2 <- "enhancers_sifoxm1_clean.txt" ## Obtained from 
sifoxm1_enhancers <- rtracklayer::import(sifoxm1_enhancers_2, format = "BED")  
sifoxm1_enhancers

##rGREAT Analysis - Obtain Exact p-values and then correct using BH / FDR:

res_senescence = great(sifoxm1_enhancers, senmayo_sasp, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_senescence = getEnrichmentTable(res_senescence) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_senescence,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment

sifoxm1_SAE <- getRegionGeneAssociations(res_senescence)
sifoxm1_SAE

all_genes <- unique(unlist(mcols(sifoxm1_SAE)$annotated_genes))

all_genes

senescence_gene_symbol <- read.table("senmayo_sasp_gene_symbol.txt", header = T, sep = "\t")
senescence_gene_symbol <- unique(senescence_gene_symbol$gene_symbol)
senescence_hits <- base::intersect(senescence_gene_symbol, all_genes)
senescence_hits


tb_senescence_df <- as.data.frame(tb_senescence)
tb_senescence_df$p_adjust <- p.adjust(tb_senescence_df$updated_pvalue, method = "BH")
write.csv(tb_senescence_df, "Senescence_H3K9me2_enhancers_sifoxm1.csv") 


res_bps = great(sifoxm1_enhancers, bps_t2g_filtered, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_bps = getEnrichmentTable(res_bps) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_bps,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment

tb_bps_df <- as.data.frame(tb_bps)
tb_bps_df$p_adjust <- p.adjust(tb_bps_df$updated_pvalue, method = "BH")
write.csv(tb_bps_df, "Inflammation_BPs_H3K9me2_Enhancers_sifoxm1.csv")


bps_gene_symbol = bps %>% dplyr::distinct(gs_name,gene_symbol) %>% as.data.frame() %>% 
  filter(gs_name %in% c("GOBP_INFLAMMATORY_RESPONSE",
                        "GOBP_CYTOKINE_PRODUCTION"))

bps_gene_symbol <- unique(bps_gene_symbol$gene_symbol)
bps_hits <- base::intersect(bps_gene_symbol, all_genes)
bps_hits


res_ap1 = great(sifoxm1_enhancers, ap1_senescence, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_ap1 = getEnrichmentTable(res_ap1) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_ap1,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment

tb_ap1_df <- as.data.frame(tb_ap1)
tb_ap1_df$p_adjust <- p.adjust(tb_ap1_df$updated_pvalue, method = "BH")
write.csv(tb_ap1_df, "AP1_Senescence_H3K9me2_Enhancers_sifoxm1.csv")


ap1_gene_symbol <- read.table("ap1_sig_gene_symbol.txt", header = T, sep = "\t")
ap1_gene_symbol <- unique(ap1_gene_symbol$gene)
ap1_hits <- base::intersect(ap1_gene_symbol, all_genes)
ap1_hits



res_hallmarks = great(sifoxm1_enhancers, hallmarks_t2g, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_hallmarks = getEnrichmentTable(res_hallmarks) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_hallmarks,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment

tb_hallmarks_df <- as.data.frame(tb_hallmarks)
tb_hallmarks_df$p_adjust <- p.adjust(tb_hallmarks_df$updated_pvalue, method = "BH")
write.csv(tb_hallmarks_df, "Hallmarks_H3K9me2_Enhancers_sifoxm1.csv")


ap1_gene_symbol <- read.table("ap1_sig_gene_symbol.txt", header = T, sep = "\t")
ap1_gene_symbol <- unique(ap1_gene_symbol$gene)
ap1_hits <- base::intersect(ap1_gene_symbol, all_genes)
ap1_hits


## Now for FOXM1 OE ------------------------------------------------------------

## Load Data for rGREAT Analysis -----------------------------------------------

## Load H3K9me2 Regions Overlapping with Enhancers obtained
## Using bedtools intersect and the ENCODE cCREs Database


dndk_enhancers_2 <- "enhancers_dndk_clean.txt"
dndk_enhancers <- rtracklayer::import(dndk_enhancers_2, format = "BED")  
dndk_enhancers

##rGREAT Analysis - Obtain Exact p-values and then correct using BH / FDR:


res_senescence = great(dndk_enhancers, senmayo_sasp, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_senescence = getEnrichmentTable(res_senescence) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_senescence,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment


tb_senescence_df <- as.data.frame(tb_senescence)
tb_senescence_df$p_adjust <- p.adjust(tb_senescence_df$updated_pvalue, method = "BH")
write.csv(tb_senescence_df, "Senescence_H3K9me2_enhancers_dndk.csv") 


res_bps = great(dndk_enhancers, bps_t2g_filtered, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_bps = getEnrichmentTable(res_bps) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_bps,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment

tb_bps_df <- as.data.frame(tb_bps)
tb_bps_df$p_adjust <- p.adjust(tb_bps_df$updated_pvalue, method = "BH")
write.csv(tb_bps_df, "Inflammation_BPs_H3K9me2_Enhancers_dndk.csv")


res_ap1 = great(dndk_enhancers, ap1_senescence, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_ap1 = getEnrichmentTable(res_ap1) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_ap1,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment

tb_ap1_df <- as.data.frame(tb_ap1)
tb_ap1_df$p_adjust <- p.adjust(tb_ap1_df$updated_pvalue, method = "BH")
write.csv(tb_ap1_df, "AP1_Senescence_H3K9me2_Enhancers_dndk.csv")

res_hallmarks = great(dndk_enhancers, hallmarks_t2g, "TxDb.Hsapiens.UCSC.hg38.knownGene")
tb_hallmarks = getEnrichmentTable(res_hallmarks) |>
  tibble::as_tibble() |>
  dplyr::select(id,genome_fraction, fold_enrichment,observed_region_hits,p_value) |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res_hallmarks,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      ))) ## This allows you to get p-values with complete accuracy, instead of 0's.
## Using BINOMIAL P-VALUE, which is appropriate for region-based enrichment

tb_hallmarks_df <- as.data.frame(tb_hallmarks)
tb_hallmarks_df$p_adjust <- p.adjust(tb_hallmarks_df$updated_pvalue, method = "BH")
write.csv(tb_hallmarks_df, "Hallmarks_H3K9me2_Enhancers_dndk.csv")
