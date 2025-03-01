library(magrittr)
library(ggrepel)
library(tidyverse)
library(vegan)
library(cluster)
library(factoextra)
library(gridExtra)
library(PerformanceAnalytics)
library(corrplot)
library(Hmisc)
library(RColorBrewer)
library(glmnet)
library(reshape2)
library(readxl)
samples_meta=read.csv("C:/My pc/Egcombio/MODA/proteomics/meta_data.csv", row.names=1, dec=",")

folder_path <-"C:/My pc/Egcombio/MODA/proteomics/GSM3375509"

# List all .txt files in the folder
file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)

# Read all files into a list of data frames
data_list <- lapply(file_list, function(file) {
  read_tsv(file, col_names = TRUE) # Read each file as a tibble
})

# Combine all data frames into one
combined_df <- bind_rows(data_list)

# Handle duplicates in EntrezGeneSymbol (e.g., by taking the mean)
combined_df <- combined_df %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(across(starts_with("M0") | starts_with("S0"), mean, na.rm = TRUE)) %>%
  ungroup()

# Filter columns that start with "S" or "M"
filtered_columns <- names(combined_df)[grepl("^S0|^M0", names(combined_df))]

# Pivot the data to wide format
wide_df <- combined_df %>%
  pivot_longer(cols = all_of(filtered_columns), names_to = "Measurement", values_to = "Value") %>%
  pivot_wider(names_from = Measurement, values_from = Value) %>%
  column_to_rownames(var = "EntrezGeneSymbol") # Set row names

data=t(wide_df)

#Remove The Zero Varaiance Proteins
varCol=apply(data, 2, var, na.rm = T)

constCol <- (varCol == 0 | is.na(varCol))

#constCol
sum(constCol)

data <- data[, !constCol]

# % of missing Values across the matrix

# Count the Sum of NAs
sum(is.na(data))
nrow(data) * ncol(data)

# The % of missingness
sum(is.na(data))/(nrow(data) * ncol(data)) * 100
round(sum(is.na(data))/(nrow(data) * ncol(data))*100,1)

# compute the missing values of Features
missingRatePerProtein = (apply(is.na(data), 2, sum)/nrow(data) ) *100

min(missingRatePerProtein)
max(missingRatePerProtein)

#Histogram for missingness rate
options(repr.plot.width=10,repr.plot.height=8)

h=hist((apply(is.na(data), 2, sum)/nrow(data) ) *100,breaks=10,
       main="Histogram for Missingness",
       xlab="percentage of missingness")
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))
abline(v = 50, col = "red", lwd = 3, lty = 2)

# checks, for each column in data, whether the percentage of missing values (NAs) is less than 50%.
good.inx=apply(is.na(data), 2, sum)/nrow(data) <0.5

data=data[,good.inx]
head(data)
dim(data)   
# Impute missing values
data.imputed=data
class(data.imputed)
head(data.imputed)
dim(data.imputed)      

range(data.imputed)    #35.4 241122.2

# Normalization:
data.imputed.logged <- log2(data.imputed + 1)
head(data.imputed.logged)
round(range(data.imputed.logged), 2)    #5.19 17.88
#Check data distrubution after and before data normalization
options(repr.plot.width=15,repr.plot.height=8)

par(mfrow=c(1,2))
plot(density( data.imputed[1,]) ,main='befor log2')
plot(density( data.imputed.logged[1,]), main='after log2' )


par(mfrow=c(1,2))
plot(density( data.imputed.logged[5,]) ,main='befor log2')
plot(density( data.imputed.logged[5,]), main='after log2' )

# Boxplot for distribution of proteins
options(repr.plot.width=10,repr.plot.height=8)
par(mar = c(8,5,2,2),mfrow=c(1,2))

# Show boxplot for the first 20 proteins before and after normalization
boxplot(data.imputed[,1:20], main="Before log2" ,horizontal=T, names=colnames(data.imputed)[1:20],las=2,col = "lightgreen")

boxplot(data.imputed.logged[,1:20], main="After log2" ,horizontal=T, names=colnames(data.imputed.logged)[1:20],las=2,col = "lightgreen")


# Boxplot for distribution of Samples
options(repr.plot.width=10,repr.plot.height=8)
par(mar = c(8,5,2,2),mfrow=c(1,2))

# Show boxplot for the first 20 Samples before and after normalizationn
boxplot(t(data.imputed[1:9,]), main="Before log2" ,horizontal=T,las=2,col = "lightgreen")

boxplot(t(data.imputed.logged[1:9,]), main="After log2" ,horizontal=T,las=2,col = "lightgreen")

#scale
data.imputed.logged.scaled=scale(data.imputed.logged,center = TRUE, scale = TRUE)
head(data.imputed.logged.scaled)

# Density plot before and after normalization & Scaling
options(repr.plot.width=20,repr.plot.height=8)

par(mar = c(8,5,2,2),mfrow=c(1,3),cex.axis=1.5)
plot(density(apply(data.imputed, 2, mean, na.rm = TRUE)),main='befor log2')
plot(density(apply(data.imputed.logged, 2, mean, na.rm = TRUE)),main='after log2')

# Boxplot for proteins Before, Normalization, Scaling
options(repr.plot.width=10,repr.plot.height=8)
par(mar = c(8,5,2,2),mfrow=c(1,3))

boxplot(data.imputed[,1:20], main="Before log2" ,horizontal=T, names=colnames(data.imputed)[1:20],las=2,col = "lightgreen")

boxplot(data.imputed.logged[,1:20], main="After log2" ,horizontal=T, names=colnames(data.imputed.logged)[1:20],las=2,col = "lightgreen")

boxplot(data.imputed.logged.scaled[,1:20], main="After log2 +scaled " ,horizontal=T,
        names=colnames(data.imputed.logged.scaled)[1:20],las=2,col = "lightgreen")

# Boxplot for proteins Samples, Normalization, Scaling

options(repr.plot.width=10,repr.plot.height=8)
par(mar = c(8,5,2,2),mfrow=c(1,3))

boxplot(t(data.imputed[1:9,]), main="Before log2" ,horizontal=T,las=2,col = "lightgreen")
boxplot(t(data.imputed.logged[1:9,]), main="After log2" ,horizontal=T,las=2,col = "lightgreen")
boxplot(t(data.imputed.logged.scaled[1:9,]), main="After log2 and scaled " ,horizontal=T,las=2,col = "lightgreen")

kmeans2 <- kmeans((data), centers = 2, nstart = 25)
# kmeans2$cluster

factoextra::fviz_cluster(kmeans2, data = data, ellipse = T,labelsize = 1)
kmeans3 <- kmeans(data, centers = 3, nstart = 25)  #DataFlair
kmeans4 <- kmeans(data, centers = 4, nstart = 25)
kmeans5 <- kmeans(data, centers = 5, nstart = 25)
#Comparing the Plots
plot1 <- fviz_cluster(kmeans2, geom = "point", data = data) + ggtitle("k = 2")
plot2 <- fviz_cluster(kmeans3, geom = "point", data = data) + ggtitle("k = 3")
plot3 <- fviz_cluster(kmeans4, geom = "point", data = data) + ggtitle("k = 4")
plot4 <- fviz_cluster(kmeans5, geom = "point", data = data) + ggtitle("k = 5")
grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

options(repr.plot.width = 20, repr.plot.height = 6)

# Ward Hierarchical Clustering
d <- dist(data, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")

plot(fit)

# split into 3 clusters

rect.hclust(fit, k=3, border="red")  # split into 13 clusters

rect.hclust(fit, k=5, border="blue")


options(repr.plot.width=8,repr.plot.height=6)
# Remove the last row from df_out to match samples_meta
df_out <- df_out[-nrow(df_out), ]
df_pca <- prcomp(data.imputed.logged,center = T,scale. = T)
df_out <- as.data.frame(df_pca$x)
# Plot PCA
library(ggplot2)

# Plot PCA
ggplot(df_out, aes(x = PC1, y = PC2, color = samples_meta$Sex, shape = samples_meta$Sex)) +
  geom_point(size = 8, alpha = 0.5) + # Size and alpha for aesthetics
  ggtitle("PCA Plot") +
  labs(color = "Sex", shape = "Sex") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 0.5, vjust = 0.5, face = "plain"),
    axis.text.y = element_text(size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),
    axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0, face = "bold"),
    axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),
    legend.title = element_blank(), # Remove legend title
    legend.text = element_text(size = 15, face = "plain"),
    strip.text = element_text(size = 15, face = "plain"),
    legend.position = "right",
    panel.background = element_rect(fill = "transparent"), # Transparent panel background
    plot.background = element_rect(fill = "transparent", color = NA), # Transparent plot background
    panel.grid.major = element_blank(), # Remove major grid
    panel.grid.minor = element_blank(), # Remove minor grid
    legend.background = element_rect(fill = "transparent"), # Transparent legend background
    legend.box.background = element_rect(fill = "transparent"), # Transparent legend box background
    axis.line = element_line(colour = "black") # Add black axis lines
  ) +
  xlab(paste0("PC 1 (", round(df_pca$sdev[1]^2 / sum(df_pca$sdev^2) * 100, 1), "%")) +
  ylab(paste0("PC 2 (", round(df_pca$sdev[2]^2 / sum(df_pca$sdev^2) * 100, 1), "%"))
Normalized_Data=t(data.imputed.logged)
library(limma)

type = as.character(samples_meta$Sex)

design <- model.matrix(~0+factor(type))

colnames(design) <- levels(factor(type))



contrast<-makeContrasts(Male-Female,levels=design)
fit <- lmFit(as.matrix(Normalized_Data), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

DEGs <- topTable(fit2,number = Inf,p.value = 0.05,adjust.method = "BH",coef = 1)
dim(DEGs)
head(DEGs)

summary(DEGs)


library(dplyr)

DEGs <- DEGs %>%
  mutate(Regulation = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not-significant"
  )
  )

head(DEGs)
library(magrittr)  # for %T>%
library(dplyr)

DEGs_up <- DEGs %>%
  filter(adj.P.Val < 0.05 & logFC > 1) %T>%
  write.csv("DEGs_UP.csv", row.names = FALSE)

DEGs_down <- DEGs %>%
  filter(adj.P.Val < 0.05 & logFC < -1) %T>%
  write.csv("DEGs_DOWN.csv", row.names = FALSE)

DEGs %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% dim %>% head %>% write.csv('limma_Proteomics.csv')
