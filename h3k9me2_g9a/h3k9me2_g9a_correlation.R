set.seed(123)

library(tidyverse)
library(ggExtra)

g9a_vs_h3k9me2 <- read.csv("g9a_vs_h3k9me2_intensity_levels.csv", header = T, row.names = NULL, sep = ";")
g9a_vs_h3k9me2 <- na.omit(g9a_vs_h3k9me2)



plot <- ggplot(g9a_vs_h3k9me2, aes(x = g9a, y = h3k9me2)) +
  geom_point(aes(shape = condition, color = condition), size = 1.5) +
  theme_bw() +
  scale_color_manual(values = c("darkgray", "black")) +
  labs(x = "G9a Intensity Levels", y = "H3K9me2 Intensity Levels") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "bottom") 
plot

plot2 <- ggMarginal(plot, type="density", groupColour = T, groupFill = T)
plot2


g9a_vs_h3k9me2_cor <- g9a_vs_h3k9me2 %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    spearman_cor = cor(g9a, h3k9me2, method = "spearman", use = "pairwise.complete.obs"),
    p_value = cor.test(g9a, h3k9me2, method = "spearman", exact = FALSE)$p.value
  )
