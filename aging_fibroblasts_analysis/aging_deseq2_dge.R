# Load Packages and Set Seed 

set.seed(123)

library(DESeq2)
library(org.Hs.eg.db)
library(ggplot2)
library(ggalt)

# Read the count matrix and metadata
counts <- read.csv("counts_matrix_young_to_old.csv", row.names = 1, sep = ";")  
# Row names are GeneIDs
coldata <- read.csv("sample_info_young_to_old.csv", row.names = 1, sep = ";")  
# Row names are Sample names

# Ensure sample names match between the two files
all(colnames(counts) == rownames(coldata))  # Should return TRUE
head(counts)



dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Condition)
dds

smallestGroupSize <- 4
keep <- rowSums(counts(dds) >=0) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "Old", "Young"),
               cooksCutoff=FALSE, independentFiltering=FALSE)
df <- as.data.frame(res)
df$gene <- row.names(df)


data = df[,"gene"]
data = as.vector(data)
annots <- AnnotationDbi::select(org.Hs.eg.db, keys=data, 
                                columns="SYMBOL", keytype="ENTREZID")
result_aging <- merge(df, annots, by.x="gene", by.y="ENTREZID")
write.csv(result_aging, "output_data\\result_aging_under40_vs_over40_cooks.csv", 
          row.names = F)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

custom_colors <- c(
  "Old" = "gray66",
  "Young" = "black")

custom_shapes <- c(
  "Old" = 17,
  "Young" = 16
)

names(custom_colors)

pca_final <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(aes(shape = Condition), size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw(base_size = 22) +
  theme(legend.position = "top",
        legend.key.size = unit(2, "point"),
        legend.key.spacing = unit(10, "pt"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values= custom_shapes) +
  geom_encircle(aes(group = Condition, color = Condition), alpha = 1, expand = 0.05, size = 1,
                s_shape = 1.2)
pca_final
saveRDS(pca_final, "pca_aging.rds")
