# Proteomics Data Analysis Pipeline

This document outlines the steps and code used for analyzing proteomics data; GSE55945(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi)
The pipeline includes data preprocessing, normalization, clustering, PCA, and differential expression analysis.

## Libraries

```R
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
```

## Data Loading and Preprocessing

### Load Metadata

```R
samples_meta <- read.csv("C:/My pc/Egcombio/MODA/proteomics/meta_data.csv", row.names = 1, dec = ",")
```

### Load and Combine Proteomics Data

```R
folder_path <- "C:/My pc/Egcombio/MODA/proteomics/GSM3375509"
file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)
data_list <- lapply(file_list, function(file) {
  read_tsv(file, col_names = TRUE)
})
combined_df <- bind_rows(data_list)
```

### Handle Duplicates and Pivot Data

```R
combined_df <- combined_df %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(across(starts_with("M0") | starts_with("S0"), mean, na.rm = TRUE)) %>%
  ungroup()

filtered_columns <- names(combined_df)[grepl("^S0|^M0", names(combined_df))]
wide_df <- combined_df %>%
  pivot_longer(cols = all_of(filtered_columns), names_to = "Measurement", values_to = "Value") %>%
  pivot_wider(names_from = Measurement, values_from = Value) %>%
  column_to_rownames(var = "EntrezGeneSymbol")
data <- t(wide_df)
```

### Remove Zero Variance Proteins

```R
varCol <- apply(data, 2, var, na.rm = TRUE)
constCol <- (varCol == 0 | is.na(varCol))
sum(constCol)
data <- data[, !constCol]
```

### Missing Value Analysis

```R
sum(is.na(data))
nrow(data) * ncol(data)
sum(is.na(data))/(nrow(data) * ncol(data)) * 100
round(sum(is.na(data))/(nrow(data) * ncol(data)) * 100, 1)

missingRatePerProtein <- (apply(is.na(data), 2, sum)/nrow(data)) * 100
min(missingRatePerProtein)
max(missingRatePerProtein)

# Histogram for missingness rate
options(repr.plot.width = 10, repr.plot.height = 8)
h <- hist((apply(is.na(data), 2, sum)/nrow(data) * 100, breaks = 10, main = "Histogram for Missingness", xlab = "Percentage of missingness")
text(h$mids, h$counts, labels = h$counts, adj = c(0.5, -0.5))
abline(v = 50, col = "red", lwd = 3, lty = 2)

# Filter out proteins with more than 50% missing values
good.inx <- apply(is.na(data), 2, sum)/nrow(data) < 0.5
data <- data[, good.inx]
head(data)
dim(data)
```

### Impute Missing Values

```R
data.imputed <- data
class(data.imputed)
head(data.imputed)
dim(data.imputed)
range(data.imputed)  # 35.4 241122.2
```

### Normalization

```R
data.imputed.logged <- log2(data.imputed + 1)
head(data.imputed.logged)
round(range(data.imputed.logged), 2)  # 5.19 17.88
```

### Check Data Distribution Before and After Normalization

```R
options(repr.plot.width = 15, repr.plot.height = 8)
par(mfrow = c(1, 2))
plot(density(data.imputed[1, ]), main = 'Before log2')
plot(density(data.imputed.logged[1, ]), main = 'After log2')

par(mfrow = c(1, 2))
plot(density(data.imputed.logged[5, ]), main = 'Before log2')
plot(density(data.imputed.logged[5, ]), main = 'After log2')
```

### Boxplot for Distribution of Proteins

```R
options(repr.plot.width = 10, repr.plot.height = 8)
par(mar = c(8, 5, 2, 2), mfrow = c(1, 2))
boxplot(data.imputed[, 1:20], main = "Before log2", horizontal = TRUE, names = colnames(data.imputed)[1:20], las = 2, col = "lightgreen")
boxplot(data.imputed.logged[, 1:20], main = "After log2", horizontal = TRUE, names = colnames(data.imputed.logged)[1:20], las = 2, col = "lightgreen")
```

### Boxplot for Distribution of Samples

```R
boxplot(t(data.imputed[1:9, ]), main = "Before log2", horizontal = TRUE, las = 2, col = "lightgreen")
boxplot(t(data.imputed.logged[1:9, ]), main = "After log2", horizontal = TRUE, las = 2, col = "lightgreen")
```

### Scaling

```R
data.imputed.logged.scaled <- scale(data.imputed.logged, center = TRUE, scale = TRUE)
head(data.imputed.logged.scaled)
```

### Density Plot Before and After Normalization & Scaling

```R
options(repr.plot.width = 20, repr.plot.height = 8)
par(mar = c(8, 5, 2, 2), mfrow = c(1, 3), cex.axis = 1.5)
plot(density(apply(data.imputed, 2, mean, na.rm = TRUE)), main = 'Before log2')
plot(density(apply(data.imputed.logged, 2, mean, na.rm = TRUE)), main = 'After log2')
```

### Boxplot for Proteins Before, After Normalization, and After Scaling

```R
options(repr.plot.width = 10, repr.plot.height = 8)
par(mar = c(8, 5, 2, 2), mfrow = c(1, 3))
boxplot(data.imputed[, 1:20], main = "Before log2", horizontal = TRUE, names = colnames(data.imputed)[1:20], las = 2, col = "lightgreen")
boxplot(data.imputed.logged[, 1:20], main = "After log2", horizontal = TRUE, names = colnames(data.imputed.logged)[1:20], las = 2, col = "lightgreen")
boxplot(data.imputed.logged.scaled[, 1:20], main = "After log2 + scaled", horizontal = TRUE, names = colnames(data.imputed.logged.scaled)[1:20], las = 2, col = "lightgreen")
```

### Boxplot for Samples Before, After Normalization, and After Scaling

```R
options(repr.plot.width = 10, repr.plot.height = 8)
par(mar = c(8, 5, 2, 2), mfrow = c(1, 3))
boxplot(t(data.imputed[1:9, ]), main = "Before log2", horizontal = TRUE, las = 2, col = "lightgreen")
boxplot(t(data.imputed.logged[1:9, ]), main = "After log2", horizontal = TRUE, las = 2, col = "lightgreen")
boxplot(t(data.imputed.logged.scaled[1:9, ]), main = "After log2 and scaled", horizontal = TRUE, las = 2, col = "lightgreen")
```

## Clustering Analysis

### K-means Clustering

```R
kmeans2 <- kmeans(data, centers = 2, nstart = 25)
kmeans3 <- kmeans(data, centers = 3, nstart = 25)
kmeans4 <- kmeans(data, centers = 4, nstart = 25)
kmeans5 <- kmeans(data, centers = 5, nstart = 25)

plot1 <- fviz_cluster(kmeans2, geom = "point", data = data) + ggtitle("k = 2")
plot2 <- fviz_cluster(kmeans3, geom = "point", data = data) + ggtitle("k = 3")
plot3 <- fviz_cluster(kmeans4, geom = "point", data = data) + ggtitle("k = 4")
plot4 <- fviz_cluster(kmeans5, geom = "point", data = data) + ggtitle("k = 5")
grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)
```

### Hierarchical Clustering

```R
d <- dist(data, method = "euclidean")
fit <- hclust(d, method = "ward")
plot(fit)
rect.hclust(fit, k = 5, border = "red")
rect.hclust(fit, k = 3, border = "red")
rect.hclust(fit, k = 5, border = "blue")
```

## Principal Component Analysis (PCA)

```R
df_pca <- prcomp(data.imputed.logged, center = TRUE, scale. = TRUE)
df_out <- as.data.frame(df_pca$x)

ggplot(df_out, aes(x = PC1, y = PC2, color = samples_meta$Sex, shape = samples_meta$Sex)) +
  geom_point(size = 8, alpha = 0.5) +
  ggtitle("PCA Plot") +
  labs(color = "Sex", shape = "Sex") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 0.5, vjust = 0.5, face = "plain"),
    axis.text.y = element_text(size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),
    axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0, face = "bold"),
    axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, face = "plain"),
    strip.text = element_text(size = 15, face = "plain"),
    legend.position = "right",
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    axis.line = element_line(colour = "black")
  ) +
  xlab(paste0("PC 1 (", round(df_pca$sdev[1]^2 / sum(df_pca$sdev^2) * 100, 1), "%")) +
  ylab(paste0("PC 2 (", round(df_pca$sdev[2]^2 / sum(df_pca$sdev^2) * 100, 1), "%"))
```

## Differential Expression Analysis

### Linear Model Fitting and Contrasts

```R
Normalized_Data <- t(data.imputed.logged)
type <- as.character(samples_meta$Sex)
design <- model.matrix(~0 + factor(type))
colnames(design) <- levels(factor(type))

contrast <- makeContrasts(Male - Female, levels = design)
fit <- lmFit(as.matrix(Normalized_Data), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

DEGs <- topTable(fit2, number = Inf, p.value = 0.05, adjust.method = "BH", coef = 1)
dim(DEGs)
head(DEGs)
summary(DEGs)
```

### Annotate DEGs

```R
DEGs <- DEGs %>%
  mutate(Regulation = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not-significant"
  ))

head(DEGs)
```

### Save DEGs

```R
DEGs_up <- DEGs %>%
  filter(adj.P.Val < 0.05 & logFC > 1) %T>%
  write.csv("DEGs_UP.csv", row.names = FALSE)

DEGs_down <- DEGs %>%
  filter(adj.P.Val < 0.05 & logFC < -1) %T>%
  write.csv("DEGs_DOWN.csv", row.names = FALSE)

DEGs %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% dim %>% head %>% write.csv('limma_Proteomics.csv')
```

---

This README file includes all the code you provided, organized into sections for clarity. Let me know if you need further adjustments! 😊
