## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(ggplot2)
library(ggimage)
library(rtracklayer)

## Import Region and Genome Annotation Data ------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene ## Genome annotation for hg38

peak_sifoxm1 <- readPeakFile("sifoxm1_peaks_signif.bed", header = F)
peak_sifoxm1
df_peak_sifoxm1 <- as.data.frame(peak_sifoxm1)

peak_foxm1_oe <- readPeakFile("dndk_peaks_signif.bed", header = F)
peak_foxm1_oe
df_peak_foxm1_oe <- as.data.frame(peak_foxm1_oe)


## FOXM1 RNAi H3K9me2 Region Characteristics -----------------------------------

peakAnno_sifoxm1 <- annotatePeak(peak_sifoxm1, tssRegion = c(-1000,100),
                                 TxDb=txdb, annoDb="org.Hs.eg.db", 
                                 overlap = "all")


peakAnno_sifoxm1_2 <- as.GRanges(peakAnno_sifoxm1)
df_peakAnno_sifoxm1 <- as.data.frame(peakAnno_sifoxm1_2)

## Several plots available:

plotAnnoPie(peakAnno_sifoxm1)
plotAnnoBar(peakAnno_sifoxm1, 
            title = "H3K9me2 Peaks Down in FOXM1 RNAi Distribution Features")
vennpie(peakAnno_sifoxm1)
plotDistToTSS(peakAnno_sifoxm1,
              title="Distribution of Lost H3K9me2 Domains\nRelative to TSS")
upsetplot(peakAnno_sifoxm1, vennpie = T) +
  theme(axis.text = element_text(size=6), 
        plot.margin = unit(c(30,30,30,30), "pt"))
upsetplot(peakAnno_sifoxm1) + theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank())
vennpie(peakAnno_sifoxm1)


## FOXM1 OE H3K9me2 Region Characteristics -----------------------------------

peakAnno_foxm1_oe <- annotatePeak(peak_foxm1_oe, tssRegion = c(-1000,100),
                              TxDb=txdb, annoDb="org.Hs.eg.db", 
                              overlap = "all")

peakAnno_foxm1_oe_2 <- as.GRanges(peakAnno_foxm1_oe)

df_peakAnno_foxm1_oe <- as.data.frame(peakAnno_foxm1_oe_2)


## Several plots available

plotAnnoPie(peakAnno_dndk)
plotAnnoBar(peakAnno_dndk, 
            title = "H3K9me2 Peaks Up in FOXM1-dNdK Distribution Features")
vennpie(peakAnno_dndk)
plotDistToTSS(peakAnno_dndk,
              title="Distribution of Gained H3K9me2 Domains\nRelative to TSS")
upsetplot(peakAnno_dndk, vennpie = T) +
  theme(axis.text = element_text(size=6), 
        plot.margin = unit(c(30,30,30,30), "pt"))
upsetplot(peakAnno_dndk) + theme(panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank())
vennpie(peakAnno_dndk)



## Overlap with cLADs and fLADs ------------------------------------------------

## Overlap with LADs, using regions available for intersection:

sifoxm1_gr <- import("sifoxm1_peaks_signif_for_overlap.bed")
lads <- import("Overlap With LADs/only_lads_hg38.bed")

# Harmonize chromosomes between query and target:

common_seqlevels <- intersect(seqlevels(sifoxm1_gr), seqlevels(lads))

sifoxm1_gr <- keepSeqlevels(sifoxm1_gr, common_seqlevels, 
                            pruning.mode="coarse")
lads <- keepSeqlevels(lads, common_seqlevels, pruning.mode="coarse")

export(sifoxm1_gr, con = "Overlap With LADs/sifoxm1_dbrs_filtered.bed", 
       format = "BED")
export(lads, con = "Overlap With LADs/lads_hg38_filtered.bed", 
       format = "BED")

sifoxm1_dbrs_filtered <- "Overlap With LADs/sifoxm1_dbrs_filtered.bed"
lads_hg38_filtered <- "Overlap With LADs/lads_hg38_filtered.bed"

# Now run enrichPeakOverlap safely:

lads_sifoxm1 <- enrichPeakOverlap(
  queryPeak     = sifoxm1_dbrs_filtered,
  targetPeak    = lads_hg38_filtered,
  TxDb          = txdb,
  pAdjustMethod = "BH",
  nShuffle      = 5000,
  chainFile     = NULL,
  verbose       = FALSE
)

write.table(lads_sifoxm1, 
            "Overlap With LADs/sifoxm1_vs_lads_overlap_statistics.txt", 
            sep = "\t")


# Harmonize chromosomes between query and target:

foxm1_oe_gr <- import("dndk_peaks_signif_for_overlap.bed")
lads <- import("Overlap With LADs/only_lads_hg38.bed")
common_seqlevels_foxm1_oe <- intersect(seqlevels(foxm1_oe_gr), seqlevels(lads))

foxm1_oe_gr <- keepSeqlevels(foxm1_oe_gr, common_seqlevels_foxm1_oe, 
                             pruning.mode="coarse")
lads_2 <- keepSeqlevels(lads, common_seqlevels_foxm1_oe, pruning.mode="coarse")

export(foxm1_oe_gr, con = "Overlap With LADs/foxm1_oe_dbrs_filtered.bed", 
       format = "BED")
export(lads_2, con = "Overlap With LADs/lads_hg38_filtered_2.bed", 
       format = "BED")

foxm1_oe_dbrs_filtered <- "Overlap With LADs/foxm1_oe_dbrs_filtered.bed"
lads_hg38_filtered_2 <- "Overlap With LADs/lads_hg38_filtered_2.bed"

# Now run enrichPeakOverlap safely:

lads_foxm1_oe <- enrichPeakOverlap(
  queryPeak     = foxm1_oe_dbrs_filtered,
  targetPeak    = lads_hg38_filtered_2,
  TxDb          = txdb,
  pAdjustMethod = "BH",
  nShuffle      = 5000,
  chainFile     = NULL,
  verbose       = FALSE
)

write.table(lads_foxm1_oe, 
            "Overlap With LADs/foxm1_oe_vs_lads_overlap_statistics.txt", 
            sep = "\t")

sessionInfo()
