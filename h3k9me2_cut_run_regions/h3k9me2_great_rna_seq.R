library(tidyverse)

set.seed(123)

## This code creates plots showing Hallmark gene sets where 
## there is enrichment in H3K9me2 changes, and plots where
## the reduction in H3K9me2 binding at enhancers are
## associated with upregulation in RNA-seq, and vice-versa:

cr_hall_sifoxm1 <- read.csv("sifoxm1/Hallmarks_H3K9me2_Enhancers_sifoxm1.csv", header=T, sep = ";")
cr_inflam_sifoxm1<- read.csv("sifoxm1/Inflammation_BPs_H3K9me2_Enhancers_sifoxm1.csv", header=T, sep = ",")
cr_sen_sifoxm1 <- read.csv("sifoxm1/Senescence_H3K9me2_enhancers_sifoxm1.csv", header=T, sep = ",")
cr_ap1_sifoxm1 <- read.csv("sifoxm1/AP1_Senescence_H3K9me2_Enhancers_sifoxm1.csv", header=T, sep = ",")

cr_inflam_dndk <- read.csv("dndk/Inflammation_BPs_H3K9me2_Enhancers_dndk.csv", header=T, sep = ",")
cr_hall_dndk <- read.csv("dndk/Hallmarks_H3K9me2_Enhancers_dndk.csv", header=T, sep = ",")
cr_sen_dndk<- read.csv("dndk/Senescence_H3K9me2_enhancers_dndk.csv", header=T, sep = ",")
cr_ap1_dndk <- read.csv("dndk/AP1_Senescence_H3K9me2_Enhancers_dndk.csv", header=T, sep = ",")


rna_seq_hall_sifoxm1 <- read.csv("sifoxm1/GSEA/siFOXM1_2024_INFLAM.csv", header=T, sep = ",")
rna_seq_inflam_sifoxm1<- read.csv("sifoxm1/GSEA/RNA SEQ siFOXM1_2024_GOBP_INFLAM.csv", header=T, sep = ",", stringsAsFactors = F)
rna_seq_sen_sifoxm1 <- read.csv("sifoxm1/GSEA/sifoxm1_sen_sasp.csv", header=T, sep = ",")
rna_seq_AP1_sifoxm1 <- read.csv("sifoxm1/GSEA/GSEA siFOXM1 AP-1 Signature.csv", header=T, sep = ",")
rna_seq_AP1_sifoxm1 <- rna_seq_AP1_sifoxm1 %>% dplyr::select(-X)

rna_seq_inflam_dndk <- read.csv("dndk/GSEA/RNA-SEQ dndk_2024_GOBP_INFLAM.csv", header=T, sep = ",")
rna_seq_hall_dndk <- read.csv("dndk/GSEA/dndk_2024_INFLAM.csv", header=T, sep = ",")
rna_seq_sen_dndk<- read.csv("dndk/GSEA/dndk_sen_sasp.csv", header=T, sep = ",")
rna_seq_AP1_dndk <- read.csv("dndk/GSEA/GSEA dNdK AP-1 Signature.csv", header=T, sep = ",")
rna_seq_AP1_dndk <- rna_seq_AP1_dndk %>% dplyr::select(-X)


cr_hall_sifoxm1 <- cr_hall_sifoxm1 %>% mutate(fold_enrichment = fold_enrichment*-1)
cr_inflam_sifoxm1<- cr_inflam_sifoxm1 %>% mutate(fold_enrichment = fold_enrichment*-1)
cr_sen_sifoxm1 <- cr_sen_sifoxm1 %>% mutate(fold_enrichment = fold_enrichment*-1)
cr_ap1_sifoxm1 <- cr_ap1_sifoxm1 %>% mutate(fold_enrichment = fold_enrichment*-1)

rna_seq_dndk <- rbind(rna_seq_inflam_dndk, rna_seq_hall_dndk, rna_seq_sen_dndk, rna_seq_AP1_dndk)
rna_seq_sifoxm1 <- rbind(rna_seq_inflam_sifoxm1, rna_seq_hall_sifoxm1, rna_seq_sen_sifoxm1, rna_seq_AP1_sifoxm1)
cr_dndk <- rbind(cr_inflam_dndk, cr_hall_dndk, cr_sen_dndk, cr_ap1_dndk)
cr_sifoxm1 <- rbind(cr_inflam_sifoxm1, cr_hall_sifoxm1, cr_sen_sifoxm1, cr_ap1_sifoxm1)

terms <- c("GOBP_CYTOKINE_PRODUCTION", "GOBP_INFLAMMATORY_RESPONSE",
           "HALLMARK_INTERFERON_ALPHA_RESPONSE",
           "HALLMARK_INTERFERON_GAMMA_RESPONSE",
           "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_INFLAMMATORY_RESPONSE",
           "HALLMARK_COMPLEMENT",
           "HALLMARK_IL2_STAT5_SIGNALING",
           "HALLMARK_IL6_JAK_STAT3_SIGNALING",
           "SASP",
           "SASP_Coppe",
           "SenMayo",
           "AP-1 Senescence Signature",
           "SASP_SenMayo", 
           "AP1_sign")

rna_seq_dndk2 <- rna_seq_dndk %>% dplyr::select(Description, NES, p.adjust) %>% 
  filter(p.adjust < 0.25) %>% 
  mutate(p.adjust = -log10(p.adjust)) %>% 
  mutate(group = "dndk RNA-seq") %>% 
  filter(Description %in% terms)


rna_seq_sifoxm1_2 <- rna_seq_sifoxm1 %>% dplyr::select(Description, NES, p.adjust) %>% 
  filter(p.adjust < 0.25) %>% 
  mutate(p.adjust = -log10(p.adjust)) %>% 
  mutate(group = "siFOXM1 RNA-seq") %>% 
  filter(Description %in% terms)


cr_dndk_2 <- cr_dndk %>% dplyr::select(id, fold_enrichment, p_adjust) %>% 
  filter(p_adjust < 0.25) %>% 
  mutate(p_adjust = -log10(p_adjust)) %>%
  filter(id %in% terms) %>% 
  mutate(group = "dndk CUT&RUN") 


cr_sifoxm1_2 <- cr_sifoxm1 %>% dplyr::select(id, fold_enrichment, p_adjust) %>% 
  filter(p_adjust < 0.25) %>% 
  mutate(p_adjust = -log10(p_adjust)) %>%
  filter(id %in% terms) %>% 
  mutate(group = "siFOXM1 CUT&RUN") 

rna_seq_dndk2 <- rna_seq_dndk2 %>% rename(Enrichment = NES, p_adjust = p.adjust)
rna_seq_sifoxm1_2 <- rna_seq_sifoxm1_2 %>% rename(Enrichment = NES, p_adjust = p.adjust)
cr_dndk_2 <- cr_dndk_2 %>% rename(Description = id, Enrichment = fold_enrichment)
cr_sifoxm1_2 <-cr_sifoxm1_2 %>% rename(Description = id, Enrichment = fold_enrichment)

all <- rbind(rna_seq_dndk2, rna_seq_sifoxm1_2, cr_dndk_2, cr_sifoxm1_2)

write.csv(all, "all_enhancers_sifoxm1_dndk.csv") 
## Added order externally and put correct names for siFOXM1 CUT&RUN

all2 <- read.csv("all_enhancers_sifoxm1_dndk.csv", header = T, sep = ",")
all2 <- dplyr::select(all2, -X)

all2 <- all2 %>% 
  dplyr::mutate(Description = case_when(
    Description == "SASP" ~ "SASP_SenMayo",
    Description == "AP1_sign" ~ "AP-1 Senescence Signature",
    .default = Description
  ))

unique(all2$Description)

order <- data.frame(
  Description = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "GOBP_INFLAMMATORY_RESPONSE",
  "GOBP_CYTOKINE_PRODUCTION",
  "SenMayo",
  "SASP_SenMayo",
  "SASP_Coppe",
  "AP-1 Senescence Signature"),
  order = 1:length(unique(all2$Description))
)


all2 <- all2 %>% 
  left_join(order, join_by(Description))

library(ggplot2)
library(forcats)

head(all2)

unique(all2$group)

all2$group <- factor(all2$group, levels = c("siFOXM1 CUT&RUN", "siFOXM1 RNA-seq", "dndk CUT&RUN", "dndk RNA-seq"))

categories <- data.frame(
  Description = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                  "HALLMARK_INFLAMMATORY_RESPONSE",
                  "HALLMARK_COMPLEMENT",
                  "HALLMARK_IL2_STAT5_SIGNALING",
                  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                  "GOBP_INFLAMMATORY_RESPONSE",
                  "GOBP_CYTOKINE_PRODUCTION",
                  "SenMayo",
                  "SASP_SenMayo",
                  "SASP_Coppe",
                  "AP-1 Senescence Signature"),
  category = c(
    "Hallmarks",
    "Hallmarks",
    "Hallmarks",
    "Hallmarks",
    "Hallmarks",
    "Hallmarks",
    "Hallmarks",
    "GO:BP",
    "GO:BP",
    "Senescence",
    "Senescence",
    "Senescence",
    "AP-1"
  )
)

all2 <- all2 %>% 
  left_join(categories, join_by(Description))

all2$category <- factor(all2$category,
                        levels = c("Hallmarks", "GO:BP", "Senescence", "AP-1"))


ggplot(all2, aes(x=group, color = Enrichment, y=fct_reorder(Description, -order))) + 
  geom_point(aes(size=p_adjust)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient2(name = "Enrichment Score", low ="deepskyblue2", mid = "white", high = "red", midpoint = 0, limits = c(-3,4.5)) +
  scale_size_continuous(name = "-log10(padjust)", range = c(2,6), limits = c(0,220), breaks = c(20,120,220)) +
  theme(axis.text = element_text(size=9), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.margin = unit(c(15,15,15,15), "pt")) +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  geom_hline(yintercept = c(1.5, 4.5, 6.5), linetype = "dashed") +
  labs(title = "CUT&RUN vs RNA-seq - Inflammation & Senescence Enhancers",
       y = "Processes", x="") +
  guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 0.4))


ggplot(all2, aes(x=group, color = Enrichment, y=fct_reorder(Description, -order))) + 
  geom_point(aes(size=p_adjust)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_color_gradient2(name = "Enrichment Score", low ="deepskyblue2", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(name = "-log10(padjust)", range = c(2,5), limits = c(0,120), breaks = c(20,60,100,140)) +
  theme(axis.text = element_text(size=9), axis.text.y = element_text(size = 8),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),plot.margin = unit(c(15,15,15,15), "pt")) +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  labs(title = "CUT&RUN vs RNA-seq - Inflammation & Senescence",
       y = "Processes", x="") +
  guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 0.4)) + facet_grid(rows = vars(category), scales = "free")


##All Hallmarks

cr__hall_dndk_2 <- cr_hall_dndk %>% dplyr::select(id, fold_enrichment, p_adjust) %>% 
  filter(p_adjust < 0.05 & fold_enrichment > 1.25) %>% 
  mutate(p_adjust = -log10(p_adjust)) %>%
  mutate(group = "dndk CUT&RUN")


cr_hall_sifoxm1_2 <- cr_hall_sifoxm1 %>% dplyr::select(id, fold_enrichment, p_adjust) %>% 
  filter(p_adjust < 0.05 & fold_enrichment < -1.25) %>% 
  mutate(p_adjust = -log10(p_adjust)) %>%
  mutate(group = "siFOXM1 CUT&RUN")


together_hallmarks <- rbind(cr_hall_sifoxm1_2, cr__hall_dndk_2)

together_hallmarks2 <- gsub("HALLMARK_"," ", together_hallmarks$id)
together_hallmarks2

together_hallmarks$id <- together_hallmarks2


together_hallmarks$group <- factor(together_hallmarks$group, levels = c("siFOXM1 CUT&RUN", "dndk CUT&RUN"))

together_hallmarks$id <- str_to_title(together_hallmarks$id)

ggplot(together_hallmarks, aes(x=group, color = fold_enrichment, y=fct_reorder(id, p_adjust))) + 
  geom_point(aes(size=p_adjust)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_color_gradient2(name = "Enrichment Score", low ="deepskyblue2", mid = "white", high = "red", midpoint = 0, limits = c(-3.5,3.5)) +
  scale_size_continuous(name = "-log10(padjust)", range = c(1,3), limits = c(0,220), breaks = c(20,120,220)) +
  theme(axis.text = element_text(size=9), axis.text.y = element_text(size = 6),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),plot.margin = unit(c(15,15,15,15), "pt")) +
  labs(title = "CUT&RUN vs RNA-seq - Inflammation & Senescence",
       y = "Processes", x="") +
  guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 0.4))


dndk_rnaseq_all_hall <- read.csv("dndk/GSEA/rna-seq - dndk_2024_HALL.csv", header = T, sep = ",") 
sifoxm1_rnaseq_all_hall <- read.csv("sifoxm1/GSEA/rna-seq - siFOXM1_2024_HALL.csv", header = T, sep = ",") 

dndk_rnaseq_all_hall2 <- dndk_rnaseq_all_hall %>% select(ID, NES, p.adjust) %>% 
  rename(id = ID, fold_enrichment = NES, p_adjust = p.adjust) %>% 
  mutate(group = "dNdK RNA-seq") %>% 
  filter(p_adjust < 0.25)

names_dndk <- gsub("HALLMARK_"," ", dndk_rnaseq_all_hall2$id)
names_dndk <- str_to_title(names_dndk)


sifoxm1_rnaseq_all_hall2 <- sifoxm1_rnaseq_all_hall %>% select(ID, NES, p.adjust) %>% 
  rename(id = ID, fold_enrichment = NES, p_adjust = p.adjust) %>% 
  mutate(group = "siFOXM1 RNA-seq") %>% 
  filter(p_adjust < 0.25)

names_sifoxm1 <- gsub("HALLMARK_"," ", sifoxm1_rnaseq_all_hall2$id)

names_sifoxm1 <- str_to_title(names_sifoxm1)


dndk_rnaseq_all_hall2$id <- names_dndk
sifoxm1_rnaseq_all_hall2$id <- names_sifoxm1

cr_names_dndk <- cr__hall_dndk_2$id
cr_names_dndk <- gsub("HALLMARK_"," ", cr_names_dndk)
cr_names_dndk <- str_to_title(cr_names_dndk)
cr_names_dndk <- as.vector(cr_names_dndk)

cr_names_sifoxm1 <- cr_hall_sifoxm1_2$id
cr_names_sifoxm1 <- gsub("HALLMARK_"," ", cr_names_sifoxm1)
cr_names_sifoxm1 <- str_to_title(cr_names_sifoxm1)

dndk_rnaseq_all_hall2 <- dndk_rnaseq_all_hall2 %>% dplyr::filter(id %in% cr_names_dndk)

sifoxm1_rnaseq_all_hall2 <- sifoxm1_rnaseq_all_hall2 %>% dplyr::filter(id %in% cr_names_sifoxm1)

all_together_rna_seq_cut_run <- rbind(together_hallmarks, dndk_rnaseq_all_hall2, sifoxm1_rnaseq_all_hall2)
all_together_rna_seq_cut_run

all_together_rna_seq_cut_run$group <- factor(all_together_rna_seq_cut_run$group, c("siFOXM1 CUT&RUN", "siFOXM1 RNA-seq",
                                                                                   "dndk CUT&RUN", "dNdK RNA-seq"))

ggplot(all_together_rna_seq_cut_run, aes(x=group, color = fold_enrichment, y=fct_reorder(id, p_adjust, .fun = max))) + 
  geom_point(aes(size=p_adjust)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient2(name = "Enrichment Score", low ="deepskyblue2", mid = "white", high = "red", midpoint = 0, limits = c(-3.5,3.5)) +
  scale_size_continuous(name = "-log10(padjust)", range = c(2,3), limits = c(0,220), breaks = c(20,120,220)) +
  theme(axis.text = element_text(size=9), axis.text.y = element_text(size = 8),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),plot.margin = unit(c(15,15,15,15), "pt")) +
  labs(title = "CUT&RUN vs RNA-seq - Hallmarks (Enhancers in H3K9me2)",
       y = "Processes", x="") +
  guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 0.4))


##New plot biased

dndk_down <- dndk_rnaseq_all_hall2 %>% filter(fold_enrichment < 0)
sifoxm1_up <- sifoxm1_rnaseq_all_hall2 %>% filter(fold_enrichment > 0)

overlap <- intersect(dndk_down$id, sifoxm1_up$id)

dndk_down <- dndk_down %>% filter(id %in% overlap)
sifoxm1_up <- sifoxm1_up %>% filter(id %in% overlap)
together_hallmarks_restrict <- together_hallmarks %>% filter(id %in% overlap)

super_restricted <- rbind(together_hallmarks_restrict, dndk_down, sifoxm1_up)
super_restricted

super_restricted$group <- factor(super_restricted$group, c("siFOXM1 CUT&RUN", "siFOXM1 RNA-seq",
                                                           "dndk CUT&RUN", "dNdK RNA-seq"))

plot1 <- ggplot(super_restricted, aes(x=group, color = fold_enrichment, y=fct_reorder(id, p_adjust))) + 
  geom_point(aes(size=p_adjust)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient2(name = "Enrichment Score", low ="deepskyblue2", mid = "white", high = "red", midpoint = 0, limits = c(-3.5,3.5)) +
  scale_size_continuous(name = "-log10(padjust)", range = c(2,4), limits = c(0,220), breaks = c(20,120,220)) +
  theme(axis.text = element_text(size=11), axis.text.y = element_text(size = 9),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),plot.margin = unit(c(15,15,15,15), "pt")) +
  labs(title = "CUT&RUN vs RNA-seq - Overlap (H3K9me2 Enhancers)",
       y = "Processes", x="") +
  guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 0.4)) +
  geom_vline(xintercept = 2.5, linetype = "dashed")

plot1 ## This plots the gene sets with opposite associations between 
      ## CUT&RUN H3K9me2 and RNA-seq - "Down H3K9me2 -> Up RNA-seq" and
      ## vice-versa

## Scatter Plot - May be Useful Later

common_hallmarks <- intersect(dndk_rnaseq_all_hall2$id, sifoxm1_rnaseq_all_hall2$id)

dndk_common <- dndk_rnaseq_all_hall2 %>% filter(id %in% common_hallmarks)
sifoxm1_common <- sifoxm1_rnaseq_all_hall2 %>% filter(id %in% common_hallmarks)

xy_plot <- common_hallmarks
xy_plot <- as.data.frame(xy_plot)
xy_plot <- xy_plot %>% rename(id = xy_plot)

sifoxm1_common <- sifoxm1_common %>% select(id, fold_enrichment) %>% 
  rename(fold_enrichment_sifoxm1 = fold_enrichment)

dndk_common <- dndk_common %>% select(id, fold_enrichment) %>% 
  rename(fold_enrichment_dndk = fold_enrichment)

xy_plot2 <- dplyr::left_join(xy_plot, sifoxm1_common, by="id")
xy_plot2 <- dplyr::left_join(xy_plot2, dndk_common, by="id")

library(ggrepel)

plot2 <- ggplot(xy_plot2, aes(x=fold_enrichment_sifoxm1, y=fold_enrichment_dndk)) +
  geom_point(shape=21, color = "blue", fill = "deepskyblue2", size = 1.5) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Significant Hallmarks - Overlap",
       x = "NES FOXM1 RNAi", y="NES FOXM1 OE") +
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "magenta", alpha = 0.1) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "magenta", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill= "green", alpha = 0.1) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill= "green", alpha = 0.1) +
  geom_label_repel(
    data=xy_plot2 %>% filter(id %in% aa), # Filter data first
    aes(label=id),
    nudge_x = 2, nudge_y = -1.5,
    label.size = 0.1, size = 3)

aaa <- xy_plot2 %>% dplyr::filter(id == " P53_pathway")

aa <- c(" Tnfa_signaling_via_nfkb", " Mtorc1_signaling", " Interferon_alpha_response", " Interferon_gamma_response", " Unfolded_protein_response")

library(patchwork)

merge <- plot2 + plot1    
merge
