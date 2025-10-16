set.seed(123)

# Load DESeq2 and other libraries
library(DESeq2)
library(org.Hs.eg.db)
library(ggplot2)
library(ggalt)

# Read the count matrix and metadata
counts <- read.csv("counts_matrix_foxm1.csv", row.names = 1, sep = ";")  # Row names are GeneIDs
coldata <- read.csv("sample_info.csv", row.names = 1, sep = ";")  # Row names are Sample names

# Ensure sample names match between the two files
colnames(counts)
row.names(coldata)
all(colnames(counts) == rownames(coldata))  # Should return TRUE
head(counts)

counts <- na.omit(counts)

str(counts)

counts_numeric <- as.data.frame(lapply(counts, function(x) as.numeric(as.character(x))))
row.names(counts_numeric) <- row.names(counts)

counts_numeric <- round(counts_numeric) ##Make sure are integer raw counts




dds <- DESeqDataSetFromMatrix(countData = counts_numeric,
                              colData = coldata,
                              design = ~ Condition) ##Start deseq2

dds
counts(dds)

smallestGroupSize <- 8
keep <- rowSums(counts(dds) >=10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))

custom_colors <- c(
  "FOXM1 OE" = "red",
  "Empty Vector" = "gray",
  "Mock" = "darkgray",
  "FOXM1 RNAi" = "deepskyblue2")

custom_shapes <- c(
  "FOXM1 OE" = 10,
  "Empty Vector" = 10,
  "Mock" = 16,
  "FOXM1 RNAi" = 16
)

names(custom_colors)

fill_data <- c("red","red", "darkgray","darkgray", "darkgray","darkgray","deepskyblue", "deepskyblue")

ages <- c("87y", "87y", "8y", "8y", "87y", "87y", "8y", "8y")

pcaData$Individual <- ages

pca_final <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(aes(shape = Condition), size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(10, "point"),
        legend.key.spacing = unit(10, "pt"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  xlim(-30,30) +
  ylim(20, -15) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values= custom_shapes) +
  geom_encircle(aes(group = Condition), alpha = 0.1, expand = 0.03, size = 1,
                s_shape = 0, fill = fill_data)
pca_final
saveRDS(pca_final, "pca_final.rds")


session_info <- devtools::session_info()
saveRDS(session_info, "session_info.RDS")
