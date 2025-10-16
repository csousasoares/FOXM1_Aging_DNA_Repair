set.seed(123)

library(factoextra)
library(FactoMineR)


##Load z-scored TPMs for each sample 

data <- read.csv("pca_2024_final_zscore.csv", row.names = 1, header = T, sep = ";")
data
head(data)


data_transposed <- t(data)
pca.data_norma <- PCA(data[,-1], scale.unit = TRUE, graph = FALSE)

pca.data <- PCA(data_transposed[-1,], scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(pca.data)
fviz_pca_ind(pca.data_norma)

var <- get_pca_var(pca.data)
var

var2 <- get_pca_var(pca.data_norma)
var2


set.seed(123)
res.km <- kmeans(var2$coord, centers = 2, nstart = 25) 

grp <- as.factor(res.km$cluster)
grp

# Sample labels for each of the 8 columns
conditions <- factor(c("mock", "mock", "sifoxm1", "sifoxm1", "vector", "vector", "dndk", "dndk"))

# Check the factor levels
levels(conditions)
# [1] "dndk" "mock" "sifoxm1" "vector"


fviz_pca_ind(pca.data, 
             col.ind = grp,
             habillage = "none",
             addEllipses = T,
             ellipse.level = 0.99,
             palette = c("Red", "Deepskyblue2"),
             ellipse.type = "confidence",
             axes = c(1,2),
             geom = c("text", "point"),
             repel = T)  +  geom_point(aes(shape = conditions, color = conditions), size=5) +
  scale_shape_manual(values = c("mock" = 16,     # solid circle
                                "sifoxm1" = 16,  # solid circle
                                "vector" = 10,   # sniper
                                "dndk" = 10)) +  # sniper
  scale_color_manual(values = c(
    "mock" = "darkgray",    
    "sifoxm1" = "deepskyblue2",  
    "vector" = "darkgray",   
    "dndk" = "red"
  )) +
  scale_x_continuous(breaks = c(-40, -30,-20,-10,0,10,20,30)) +
  scale_y_continuous(breaks = c(-40, -30,-20,-10,0,10,20,30)) +
  theme()


session_info <- devtools::session_info()
saveRDS(session_info, "session_info.RDS")
