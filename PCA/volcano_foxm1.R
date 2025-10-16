library(tidyverse)
library(ggplot2)
library(ggrastr)

deseq2_table <- read.csv("dndk_volcano.csv", header = T, sep = ";") ##dndk or sifoxm1


deseq2_table <- na.omit(deseq2_table)

deseq2_table$padj <- -log10(deseq2_table$fdr)

deseq2_table_filtered_2 <- deseq2_table %>% 
  mutate(direction = case_when(
  log2fc > 0 & fdr < 0.05 ~ "Upregulated",
  log2fc < 0 & fdr < 0.05 ~ "Downregulated",
  .default = "No Change"
))

colors <- c(
  "Downregulated" = "deepskyblue2",
  "No Change" = "gray",
  "Upregulated" = "red"
)

ggplot(deseq2_table_filtered_2, aes(x = log2fc, y = padj, colour = direction), alpha = 0.8) +
  rasterise(geom_point(shape = 16, size = 0.8), dpi = 300) +
  theme_bw() +
  xlim(-7, 7) +
  scale_color_manual(values = colors, name = "Direction") +
  geom_vline(xintercept = c(0), linetype = "dashed") +
  geom_hline(yintercept = 1.301, linetype = "dashed") +
  labs(y = "-log10(padjust)", x = "Log2FC") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

length(which(deseq2_table_filtered_2$direction=="Downregulated"))  
length(which(deseq2_table_filtered_2$direction=="Upregulated"))  
length(which(deseq2_table_filtered_2$direction=="No Change"))  


deseq2_table <- read.csv("sifoxm1_volcano.csv", header = T, sep = ";") ##dndk or sifoxm1


deseq2_table <- na.omit(deseq2_table)

deseq2_table$padj <- -log10(deseq2_table$fdr)

deseq2_table_filtered_2 <- deseq2_table %>% 
  mutate(direction = case_when(
    log2fc > 0 & fdr < 0.05 ~ "Upregulated",
    log2fc < 0 & fdr < 0.05 ~ "Downregulated",
    .default = "No Change"
  ))

colors <- c(
  "Downregulated" = "deepskyblue2",
  "No Change" = "gray",
  "Upregulated" = "red"
)

ggplot(deseq2_table_filtered_2, aes(x = log2fc, y = padj, colour = direction), alpha = 0.8) +
  rasterise(geom_point(shape = 16, size = 0.8), dpi = 300) +
  theme_bw() +
  xlim(-4, 5) +
  scale_color_manual(values = colors, name = "Direction") +
  geom_vline(xintercept = c(0), linetype = "dashed") +
  geom_hline(yintercept = 1.301, linetype = "dashed") +
  labs(y = "-log10(padjust)", x = "Log2FC") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

length(which(deseq2_table_filtered_2$direction=="Downregulated"))  
length(which(deseq2_table_filtered_2$direction=="Upregulated"))  
length(which(deseq2_table_filtered_2$direction=="No Change"))  

