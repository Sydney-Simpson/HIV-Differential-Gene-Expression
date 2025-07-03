Data Analysis
================
2025-06-17

``` r
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warnings = FALSE)
```

``` r
library(dplyr)
library(gprofiler2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(sqldf)
library(limma)
library(dplyr)
library(tidyr)
library(gprofiler2)
library(PCAtools)
library(ggplotify)
library(gridExtra)
library(grid)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(org.Hs.eg.db)

dat <- read.table("Processed_data.txt",  sep = "\t")
sample_info <- read.table("all_sample_info.txt",  sep = "\t")
```

``` r
# Filter for samples from negative HIV patients
negatives <- subset(sample_info, status == "HIV negative")

# Pick final samples from patients before they tested positive (for GSE195434 data)
selected_negs <- negatives %>%
  arrange(patient_id, time) %>%  # Order by patient id and time
  group_by(patient_id) %>%       # Group by patient id
  slice_max(order_by = time, n=1) %>%      #Select the last sample
  ungroup()                             # Remove Grouping

# Filter for samples from positive patients
positives <- subset(sample_info, status %in% c( "HIV-1 infected", "HIV-1 and HIV-2 dual infected", "HIV-2 infected"))

# Pick first samples from patients when they tested positive- and were treatment naive
selected_pos <- positives %>%
  arrange(patient_id, time) %>%    # Order by patient id and time
  group_by(patient_id) %>%        # Group by patient id
  slice_min(order_by = time, n =1) %>%  #Select the earliest sample
  ungroup()                             # Remove Grouping

# Combine negative and positive samples
selected_sample_info <-rbind(selected_negs, selected_pos) 

# Extract GEO accession numbers for selected samples
selected_samples <- selected_sample_info$geo_accession

# Subset the data to include only the selected samples
selected_data <- dat[,c(selected_samples)]
```

``` r
# Select only the last pre infection and the first post infection sample taken from each patient in GSE195434
selected_negs <- sqldf("WITH negatives AS (
                      SELECT * 
                      FROM sample_info
                      WHERE status = 'HIV negative'
                    ),
                   ranked_samples AS( SELECT *,
                     ROW_NUMBER() OVER (PARTITION BY patient_id ORDER BY time DESC) as rn
                     FROM negatives
                   )
                   SELECT *
                   FROM ranked_samples
                   WHERE rn=1
                   ")

# Pick first samples from patients when they tested positive- and were treatment naive
selected_pos <- sqldf("WITH positives AS (
                      SELECT * 
                      FROM sample_info
                      WHERE NOT status ='HIV negative'
                    ),
                   ranked_samples AS( SELECT *,
                     ROW_NUMBER() OVER (PARTITION BY patient_id ORDER BY time ASC) as rn
                     FROM positives
                   )
                   SELECT *
                   FROM ranked_samples
                   WHERE rn=1
                   ")


# Combine selected negative and positive samples
selected_sample_info <- sqldf("SELECT * FROM selected_negs 
                              UNION
                              SELECT * FROM selected_pos
                              ")

# Extract GEO accession numbers for selected samples
selected_samples <- sqldf("SELECT geo_accession 
                          FROM selected_sample_info")

selected_samples <- selected_samples$geo_accession

# SQL does not support dynamic column selection so we use base R to index our columns of interest
selected_data <- dat[,c(selected_samples)]
```

``` r
boxplot(selected_data)
```

![](Analysis_files/figure-gfm/get%20a%20quick%20look%20at%20the%20distribution-1.png)<!-- -->

``` r
# Select the status column from the selected sample information
labels <- dplyr::select(selected_sample_info, status)

# Transpose the selected data and convert it to a dataframe
transposed_data <- data.frame(
  t(
    as.matrix(selected_data)
    )
  )

# Label the transposed data
labelled_data <- cbind(labels, transposed_data)

# Calculate mean values for each group
mean_values <- labelled_data %>%
  group_by(status) %>%
  dplyr::summarize(across(everything(), \(x) mean(x, na.rm = TRUE))) %>% # Compute mean for each column, remove any NAs
  tibble::column_to_rownames(var = "status") # Convert the status column to row names
```

``` r
# Define groups based on disease status
groups <- factor(
  unique(labelled_data$status))

# Create a design matrix for linear modelling
design <- model.matrix(~status + 0, labelled_data)

# Clean column names and remove the 'status' prefix
colnames(design) <- make.names(colnames(design))
colnames(design) <- gsub("status", "", colnames(design))
colnames(design) <- factor(colnames(design), levels = colnames(design))

# Transpose the data for linear modeling and fit the model
t_labelled_data <- t(labelled_data[,-1])
fit <- lmFit(t_labelled_data, design)

# Set up contrasts of interest for comparison
coi <- c("HIV.1.infected-HIV.negative", "HIV.2.infected-HIV.negative", "HIV.1.and.HIV.2.dual.infected-HIV.negative")
contrasts <- makeContrasts(contrasts =coi,
                           levels = design)

# Recalculate model coefficients based on contrasts
fit2 <- contrasts.fit(fit, contrasts)

# Apply empirical Bayes statistics to the model
fit2 <- eBayes(fit2)

# Summarize and visualize results
results <- decideTests(fit2)

summary(results)
```

    ##        HIV.1.infected-HIV.negative HIV.2.infected-HIV.negative
    ## Down                             4                           0
    ## NotSig                       21066                       21186
    ## Up                             116                           0
    ##        HIV.1.and.HIV.2.dual.infected-HIV.negative
    ## Down                                            0
    ## NotSig                                      21185
    ## Up                                              1

``` r
vennDiagram(results,
            cex = c(0.8, 1.0, 1.0),
            names = c("HIV-1 vs HIV negative", "HIV-2 vs HIV negative", "Dual infected vs HIV negative"))
```

![](Analysis_files/figure-gfm/Analyze%20different%20groups%20using%20limma-1.png)<!-- -->

``` r
# Get the adjusted p value for each coefficient
tT_HIV1 <- topTable(fit2, coef =1, adjust="BH", n = Inf)
tT_HIV2 <- topTable(fit2, coef =2, adjust="BH", n = Inf)
tT_dual <- topTable(fit2, coef =3, adjust="BH", n = Inf)

# Extract adjusted p values for each contrast
adj_p <- list()

for (contrast in coi) {
  topTable_results <- topTable(fit2,
                               coef = contrast,
                               n= Inf,
                               adjust = "BH")
  
  adj_p[[contrast]] <- data.frame(row.names = rownames(topTable_results), 
                                  adj.P.Val = topTable_results$adj.P.Val)
}
```

``` r
# Transpose data frame of mean values
t_means <- data.frame(t(mean_values))

# Calculate the average fold change for HIV positive and negative samples for each gene
change <- t_means - t_means$HIV.negative

# Remove the HIV.negative column, as it is now zero
change <- data.frame(change[,-1]) 

# Merge adjusted p-values for all contrasts
adjusted_ps <- merge(adj_p[[1]], adj_p[[2]], by = "row.names") %>%
  tibble::column_to_rownames(var = "Row.names")

# Continue merging with the third set of p-values
adjusted_p_values <- data.frame(merge(adjusted_ps, adj_p[[3]], by = "row.names")) %>%
  tibble::column_to_rownames(var = "Row.names")

# Clean column names by removing '-HIV.negative'
colnames(adjusted_p_values) <- gsub("-HIV.negative", "",names(adj_p))

#Save the adjusted p-values table to a file
write.table(adjusted_p_values, "adjusted_p_values.txt", row.names= TRUE, quote= FALSE, sep="\t")

# Create data frame with fold change and adjusted p-value data for each group of HIV infected samples

# Check if row names are identical, if they are combine the columns
if (identical(rownames(change), rownames(adjusted_p_values))) {
  dual <- data.frame(Change = change$HIV.1.and.HIV.2.dual.infected, pvalue = adjusted_p_values$HIV.1.and.HIV.2.dual.infected, 
                     row.names = rownames(change))
  HIV1 <- data.frame(Change = change$HIV.1.infected, pvalue = adjusted_p_values$HIV.1.infected, 
                     row.names = rownames(change))
  HIV2 <- data.frame(Change = change$HIV.2.infected, pvalue = adjusted_p_values$HIV.2.infected, 
                     row.names = rownames(change))
} else {
  # If they don't match then order them accordingly
  change <- change[order(row.names(adjusted_p_values)),]
  
  dual <- data.frame(Change = change$HIV.1.and.HIV.2.dual.infected, pvalue = adjusted_p_values$HIV.1.and.HIV.2.dual.infected, 
                     row.names = rownames(change))
  HIV1 <- data.frame(Change = change$HIV.1.infected, pvalue = adjusted_p_values$HIV.1.infected, 
                     row.names = rownames(change))
  HIV2 <- data.frame(Change = change$HIV.2.infected, pvalue = adjusted_p_values$HIV.2.infected, 
                     row.names = rownames(change))
}

# Store datasets in a list
datasets <- list()
datasets[["HIV1"]] <- HIV1
datasets[["HIV2"]] <- HIV2
datasets[["Dual Infected"]] <- dual

# Plot each condition for a quick examination
for (i in 1:length(datasets)){
  dat <- datasets[[i]] 
  dat <- dat %>% 
    mutate(Significant = pvalue < 0.05) 
  
  # Add labels genes with the lowest 20 p values, and most significant changes in expression
  p_labels <- rownames(dat[order(dat$pvalue),])[1:20]
  fold_labels <- rownames(dat[order(-abs(dat$Change)),])[1:20]
  
  dat$gene <- rownames(dat)
  dat <- dat %>%
    mutate(PLabel = ifelse(gene %in% p_labels, gene, "")) %>%
    mutate(CLabel = ifelse(gene %in% fold_labels, gene, ""))
  
  plot1 <- ggplot(data = dat, mapping = aes(x = Change, y = -log10(pvalue), label = PLabel, col =  Significant)) +
    geom_point() +
    labs(title = paste(names(datasets)[i], "Genes with Statistically Significant Differential Expression" ), y = "-log10(adjusted p value)", x = "log2 difference") +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_text_repel(col="black")
  
  plot2 <- ggplot(data = dat, mapping = aes(x = Change, y = -log10(pvalue), label = CLabel, col =  Significant)) +
    geom_point() +
    labs(title = paste(names(datasets)[i], "Most Differentially Expressed Genes" ), y = "-log10(adjusted p value)", x = "log2 difference") +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_text_repel(col="black")
    
  print(plot1)
  print(plot2)
}
```

    ## Warning: ggrepel: 6 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-1.png)<!-- -->

    ## Warning: ggrepel: 6 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-2.png)<!-- -->

    ## Warning: ggrepel: 20 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-3.png)<!-- -->

    ## Warning: ggrepel: 17 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-4.png)<!-- -->

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-5.png)<!-- -->

    ## Warning: ggrepel: 17 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-6.png)<!-- -->

``` r
# The HIV-2 Infected and Dual Infected conditions only have 1 statistically significant gene in total likely due to low sample numbers. However, the HIV-1 infected samples are very interesting so we'll save those graphs

# Add labels genes with the lowest 20 p values, and most significant changes in expression
p_labels <- rownames(HIV1[order(HIV1$pvalue),])[1:20]
fold_labels <- rownames(HIV1[order(-abs(HIV1$Change)),])[1:20]
  
HIV1$gene <- rownames(HIV1)
HIV1 <- HIV1 %>%
  mutate(Significant = pvalue < 0.05) %>%
  mutate(PLabel = ifelse(gene %in% p_labels, gene, "")) %>%
  mutate(CLabel = ifelse(gene %in% fold_labels, gene, ""))

plot1 <- ggplot(data = HIV1, mapping = aes(x = Change, y = -log10(pvalue), label = PLabel, col =  Significant)) +
  geom_point() +
  labs(title = "Genes with the most statistically significant differential \n expression between HIV-1 infected and HIV negative individuals", y = "-log10(adjusted p value)", x = "log2 difference", color = "Statistically Significant (p < 0.05)")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(col="black")
  
plot2 <- ggplot(data = HIV1, mapping = aes(x = Change, y = -log10(pvalue), label = CLabel, col =  Significant)) +
  geom_point() +
  labs(title = "Genes with the highest fold change between \n HIV-1 infected and HIV negative individuals", y = "-log10(adjusted p value)", x = "log2 difference", color = "Statistically Significant (p < 0.05)") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(col="black")
    
print(plot1)
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-7.png)<!-- -->

``` r
print(plot2)
```

    ## Warning: ggrepel: 11 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-8.png)<!-- -->

``` r
# There was only one gene that had statistically significant differential expression for dual infection
dual$gene <- rownames(dual)
dual <- dual %>%
  mutate(Significant = pvalue < 0.05) %>%
  mutate(PLabel = ifelse(gene %in% c("TMEM119"), gene, ""))

# Show the plot
dual_infect_p <- ggplot(data = dual, mapping = aes(x = Change, y = -log10(pvalue), label = PLabel, col =  Significant)) +
  geom_point() +
  labs(title = "Differential expression between HIV negative /n and dual infected individuals", y = "-log10(adjusted p value)", x = "log2 difference", color = "Statistically Significant (p < 0.05)")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(col="black")

# Make a new data frame for this graph, and substitute a new value for dual infected individuals
TMEM_data <- dplyr::select(labelled_data, status, TMEM119) 
TMEM_data$status <- gsub("HIV-1 and HIV-2 dual infected", "Dual infected", TMEM_data$status) 
TMEM_data$status <- factor(TMEM_data$status, levels = c(
    "HIV negative", "HIV-1 infected", "Dual infected", "HIV-2 infected"))

TMEM_plot <- ggplot(TMEM_data ,aes(x = status, y = TMEM119, group = status, fill = status)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize = 1 ) +
  geom_boxplot()+
  labs(y = "TMEM119 log2 expression", x = NULL)

TMEM_plot
```

![](Analysis_files/figure-gfm/Identify%20changes%20in%20gene%20expression%20then%20match%20with%20p%20values%20for%20each%20infection%20status-9.png)<!-- -->

``` r
# Filter genes based on statistical significance
sig_genes <- rownames(subset(adjusted_p_values,
              HIV.1.infected <= 0.05 | HIV.2.infected <= 0.05 | HIV.1.and.HIV.2.dual.infected <= 0.05))

# Filter genes based on change over 1.5x. Note that this data is log2 transformed and log2(1.5) is 0.6
change2 <- change
changed_genes <- row.names(subset(change2,
              abs(HIV.1.infected) >= 0.6 | abs(HIV.2.infected) >= 0.6 | abs(HIV.1.and.HIV.2.dual.infected) >= 0.6)) 

# Get a list of unique genes that met either statistical significance or gold change criteria
unique_genes <- unique(c(sig_genes, changed_genes))

# Ensure the rownames, aka GeneSymbols in your selected data is in the same format as the list of unique genes of interest
rownames(selected_data) <- make.names(rownames(selected_data))

# Select genes of interest from selected sample data
selected_genes <- selected_data[unique_genes,]

# Look at distribution of data in selected dataset
boxplot(selected_genes)
```

![](Analysis_files/figure-gfm/create%20heatmaps%20and%20narrow%20down%20genes%20list-1.png)<!-- -->

``` r
# Get relevant information for each sample
sample_groups <- dplyr::select(selected_sample_info, geo_accession, status, dataset) %>%
  tibble::column_to_rownames(var = "geo_accession")

# Remove outliers------------------------------------------------------------------------
# Create a function to find outliars
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# Start a list for outliers
outliers <- list()

# Create a table for this analysis
table <- t(selected_genes)
table <- merge(sample_groups, table, by = "row.names") 
table <- table %>%
  tibble::column_to_rownames(var = "Row.names")
columns <- names(table)[sapply(table, is.numeric)]


# Loop through each quantitative column and make a list of outliers for each group
for (col in columns) {
  # Select relevant data
  dat <- dplyr::select(table, col, status)%>%
    mutate(sample = row.names(table))
  # Make a dummy column, since mutate can't use variable names for columns
  dat$value <- dat[,1]
  # Identify outliers
  dat <- dat %>%
        group_by(status) %>%
        mutate(outlier = case_when(findoutlier(value) ~ sample))
  #outlying_samples <- dat$outlier[!is.na(dat$outlier)]
  outliers <- c(outliers, dat$outlier[!is.na(dat$outlier)])
}
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(col)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(col))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# Create a table to view how many times each sample was an outlier
outlier_table <- rev(stack(table(unlist(outliers)))) %>%
  arrange(desc(values))

print(head(outlier_table))
```

    ##          ind values
    ## 1 GSM6038882     76
    ## 2 GSM6038877     42
    ## 3 GSM5836724     32
    ## 4 GSM5836696     22
    ## 5 GSM5836653     16
    ## 6 GSM6038895     16

``` r
# Select samples to remove, here we are removing any where 40 or more genes out of our 225 genes of interest were outliers
samples_to_rm <- droplevels(outlier_table$ind[which(outlier_table$values >= 40)])

# Show which outliers are being removed
col_vector <- rep("lightgray", length(colnames(selected_genes)))
col_vector[which(colnames(selected_genes) %in% samples_to_rm)] <- "red"
  
boxplot(selected_genes, col = col_vector)
```

![](Analysis_files/figure-gfm/create%20heatmaps%20and%20narrow%20down%20genes%20list-2.png)<!-- -->

``` r
selected_genes <- selected_genes[,!(colnames(selected_genes) %in% samples_to_rm)]
sample_groups <- sample_groups[!(row.names(sample_groups) %in% samples_to_rm),]
```

``` r
# Specify colors for heatmap
my_colors <- list(
  status = c("HIV-1 infected" = "red", "HIV negative" = "green4", "HIV-1 and HIV-2 dual infected" = "purple", "HIV-2 infected" = "blue"),
             dataset = c(GSE195434 = "black", GSE200606 = "white")
  )

# First, look at samples and see if they are grouped by dataset or condition
corMat <- cor(selected_genes, use = "c")
corMat <- pheatmap(corMat,
                   annotation_col = sample_groups,
                   annotation_colors = my_colors,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   main = "Correlation Between Samples")

print(corMat)
```

![](Analysis_files/figure-gfm/create%20heat%20maps-1.png)<!-- -->

``` r
# Now look at genes of interest, remove the dataset information from the annotation
sample_groups <- sample_groups[, !(names(sample_groups) %in%  c("dataset")), drop = FALSE]

# Convert data values to show difference from the mean
selected_genes_diff <- selected_genes - rowMeans(selected_genes)

all_genes <- pheatmap(selected_genes_diff, 
                      annotation_col = sample_groups,
                      annotation_colors = my_colors,
                      show_rownames = FALSE,
                      show_colnames = FALSE,
                      main = "Heatmap of Selected Genes")
print(all_genes)
```

![](Analysis_files/figure-gfm/create%20heat%20maps-2.png)<!-- -->

``` r
# What if we limit this to the top 10 genes?
top10 <- topTable(fit2, adjust = "BH")
top10_data <- selected_genes_diff[rownames(top10),]


top10_p <- pheatmap(top10_data, 
                    annotation_col = sample_groups,
                    annotation_colors = my_colors,
                    show_colnames = FALSE,
                    main = "10 most statistically significant genes")

print(top10_p)
```

![](Analysis_files/figure-gfm/create%20heat%20maps-3.png)<!-- -->

``` r
# What if we limit this to the top 10 genes in regard to fold change?
# Make a new dataframe that shows the max absolute fold change
change$gene <- row.names(change)
change <- change %>%
  rowwise()%>%
  mutate(max = pmax(HIV.1.and.HIV.2.dual.infected, HIV.1.infected, HIV.2.infected, na.rm = TRUE))%>%
  mutate(min = pmin(HIV.1.and.HIV.2.dual.infected, HIV.1.infected, HIV.2.infected, na.rm = TRUE))%>%
  mutate(FC = case_when(abs(min) > abs(max) ~ min, 
                        abs(min) < abs(max) ~ max))

change <- change[order(-abs(change$FC)),]
row.names(change) <- change$gene
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
top10_changed_genes <- change[1:10,]
top10_diff <- selected_genes_diff[top10_changed_genes$gene,]

top10_diff <- pheatmap(top10_diff, 
                       annotation_col = sample_groups,
                       annotation_colors = my_colors,
                       show_colnames = FALSE,
                       main = "Top 10 genes with the highest fold changes")

print(top10_diff)
```

![](Analysis_files/figure-gfm/create%20heat%20maps-4.png)<!-- -->

``` r
# Select genes with a fold change > 1 and a statistically significant p value
genes_list <- list()
for (i in 1:length(datasets)){
  dat <- datasets[[i]]
  genes_to_add <- rownames(subset(dat,
              pvalue <= 0.05 & abs(Change) >= 1))
  genes_list <- c(genes_list, genes_to_add)
}

Sig_Change <- selected_genes_diff[unlist(genes_list),]
most_sig <- pheatmap(Sig_Change,
                     annotation_color = my_colors,
                     annotation_col = sample_groups,
                     show_colnames = FALSE,
                     main = "Genes with statistical signifcance and \n at least 1 fold change")

print(most_sig)
```

![](Analysis_files/figure-gfm/create%20heat%20maps-5.png)<!-- -->

``` r
# Transpose data
PCA_data <- selected_genes[unlist(genes_list),]
  
pca_sample_groups <- sample_groups[colnames(selected_genes),, drop = FALSE]

# Add in a column to preserve the type of HIV infection, then simplify status to show if someone is positive or not
pca_sample_groups$type <- pca_sample_groups$status
pca_sample_groups$type <- gsub("HIV-1 infected", "HIV-1",
                            gsub("HIV-1 and HIV-2 dual infected", "HIV-1 and HIV-2",
                            gsub("HIV-2 infected","HIV-2",  pca_sample_groups$type)))

pca_sample_groups$status <- gsub("HIV-1 infected", "HIV positive",
                            gsub("HIV-1 and HIV-2 dual infected", "HIV positive",
                            gsub("HIV-2 infected","HIV positive",  pca_sample_groups$status)))

pca <- pca(PCA_data, metadata = pca_sample_groups, center = TRUE, scale = TRUE)

# Make a biplot with ellipses showing sample populations
biplot(pca, 
       colby = "status", colkey = c("HIV positive" = "#FC4E07", "HIV negative" ="#00AFBB"),
       shape = "type",
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 0,
       title = "Principal Component Analysis HIV infected vs HIV negative",
       legendPosition = 'right', legendLabSize = 10, legendIconSize = 3.0,
       lab = ""
       )
```

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-1.png)<!-- -->

``` r
# Make a biplot that shows important loading info
biplot(pca, 
       colby = "status", colkey = c("HIV positive" = "#FC4E07", "HIV negative" = "#00AFBB"),
       showLoadings = TRUE,
       sizeLoadingsNames = 4,
       legendPosition = 'right', legendLabSize = 10, legendIconSize = 3.0,
       lab = ""
       )
```

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_segment()`).

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_label_repel()`).

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-2.png)<!-- -->

``` r
# Use the elbow method to identify which PCs to keep
elbow <- findElbowPoint(pca$variance)

# Make a screeplot showing results of the elbow method, as well as which components explain 90% of variance
screeplot(pca,
          components = getComponents(pca, 1:20),
          title = "Screeplot",
          hline = 90,
          vline = elbow)+
  annotate( "text",x = elbow + 2.5, y = 25,
                 label = "Elbow method", vjust = -1, size = 5)+
  annotate("text", x = 9, y = 80,
  label = "Explains 90% of variance", size = 5
  )
```

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-3.png)<!-- -->

``` r
# plot loadings
plotloadings(pca, caption = "Top 5% of variables")
```

    ## Warning: ggrepel: 14 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-4.png)<!-- -->

``` r
plotloadings(pca,
  components = getComponents(pca, c(1:3)),
  drawConnectors = TRUE, labSize = 3, caption = "Top 5% of variables")
```

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-5.png)<!-- -->

``` r
plotloadings(pca,
  components = getComponents(pca, c(1)),
  drawConnectors = TRUE, labSize = 3, caption = "Top 5% of variables",
  legendPosition = 'none') + coord_flip()
```

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-6.png)<!-- -->

``` r
plotloadings(pca,
  components = getComponents(pca, c(2)),
  drawConnectors = TRUE, labSize = 3, caption = "Top 5% of variables") + coord_flip()
```

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-7.png)<!-- -->

``` r
plotloadings(pca,
  components = getComponents(pca, c(3)),
  drawConnectors = TRUE, labSize = 3, caption = "Top 5% of variables") + coord_flip()
```

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-8.png)<!-- -->

``` r
# Create a pairs plot
pairsplot(pca,
          components = getComponents(pca, 1:3),
          colby= "status",
          trianglelabSize = 10,
          hline = 0, vline = 0,
          triangle =  TRUE,
          gridlines.major = FALSE,
          gridlines.minor = FALSE)
```

![](Analysis_files/figure-gfm/PCA%20for%20genes%20with%20p%20value%20less%20than%20or%20equal%20to%200.05%20and%20a%201%20or%20greater%20fold%20change-9.png)<!-- -->

``` r
# Create a list the genes isolated from PCA
top <- list("EPSTI1", "IFI44", "IFI44L", "LY6E", "XAF1", "CD52", "EOMES", "GZMH", "IFI27", "HERC5", "IFIT1", "RSAD2", "NKG7")

# Select these genes from the data frame showing the differences in expression for each sample to the mean
top_data <- selected_genes_diff[unlist(top),]

# Create a heatmap
top_genes <- pheatmap(top_data,
                     annotation_color = my_colors,
                     annotation_col = sample_groups,
                     show_colnames = FALSE,
                     main = "Top Genes from PCA")

top_genes
```

![](Analysis_files/figure-gfm/select%20top%20genes%20from%20PCA%20and%20create%20heatmap-1.png)<!-- -->

``` r
# Select these genes from the expression matrix, remember to remove the samples denoted as outliers
top_exprs <- labelled_data[!(row.names(labelled_data) %in% samples_to_rm),c(unlist(top), "status")]
# Assuming your data frame is named `df`, and the qualitative column is named `Qualitative`
columns <- names(top_exprs)[sapply(top_exprs, is.numeric)]

# Write functions to tell if a value is an outlier
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# Clean up labels
top_exprs$status <- gsub("HIV-1 and HIV-2 dual infected", "Dual infected", top_exprs$status) 
top_exprs$status <- factor(top_exprs$status, levels = c(
    "HIV negative", "HIV-1 infected", "Dual infected", "HIV-2 infected"))

# Start a list for outliers
outliers <- list()

# Loop through each quantitative column and create a plot
for (col in columns) {
  # Select relevant data
  dat <- dplyr::select(top_exprs, col, status)%>%
    mutate(sample = row.names(top_exprs))
  # Make a dummy column, since mutate can't use variable names for columns
  dat$value <- dat[,1]
  # Identify outliers
  dat <- dat %>%
        group_by(status) %>%
        mutate(outlier = case_when(findoutlier(value) ~ sample))
  #outlying_samples <- dat$outlier[!is.na(dat$outlier)]
  outliers <- c(outliers, dat$outlier[!is.na(dat$outlier)])
  # Create plot
  p <- ggplot(dat, aes_string(x = 'status', y = col, color = 'status')) +
    geom_boxplot() +
    geom_point()+
    labs(title = paste(col, "Expression by Infection Status")) +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3, size = 3)+
    scale_color_manual(values = c("green4", "red", "purple", "blue"))
  
  print(p)
}
```

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-1.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-2.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-3.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-4.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-5.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-6.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-7.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-8.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-9.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-10.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-11.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-12.png)<!-- -->![](Analysis_files/figure-gfm/create%20graphs%20to%20show%20distribution%20for%20each%20gene%20of%20interest-13.png)<!-- -->

``` r
# How many times is the same sample and outlier in the dataset?
outlier_table <- rev(stack(table(unlist(outliers)))) %>%
  arrange(desc(values))

print(head(outlier_table))
```

    ##          ind values
    ## 1 GSM5836663      5
    ## 2 GSM5836707      4
    ## 3 GSM5836724      4
    ## 4 GSM6038904      4
    ## 5 GSM5836696      3
    ## 6 GSM6038896      3

``` r
# Select the list of unique genes from the HIV-1 dataframe which contained log2 fold change values
unique_HIV1 <- HIV1[unique_genes,]%>%
  mutate(gene = unique_genes)
gene_list <- unique_HIV1$Change
names(gene_list) <- unique_HIV1$gene

#sort by fold change
gene_list <- sort(gene_list, decreasing = TRUE)

# Get info
gse <- gseGO(geneList = gene_list,
             ont ="ALL",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH")

# Make plots
enrichplot::dotplot(gse, split= ".sign", showCategory = 8) + facet_grid(.~.sign)
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment-1.png)<!-- -->

``` r
gse2 <- pairwise_termsim(gse)

emapplot(gse2, showCategory = 10)
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment-2.png)<!-- -->

``` r
cnetplot(gse, categorySize="pvalue", showCategory = 5, foldChange=gene_list) 
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment-3.png)<!-- -->

``` r
heatplot(gse, foldChange = gene_list, showCategory = 10, ) + scale_fill_gradient2(low = "red", , high = "blue", midpoint = 0)
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment-4.png)<!-- -->

``` r
# Make plots
enrichplot::dotplot(gse, split= ".sign", showCategory = 8) + facet_grid(.~.sign)
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment%20smaller%20plots%20for%20larger%20text-1.png)<!-- -->

``` r
gse2 <- pairwise_termsim(gse)

emapplot(gse2, showCategory = 10)
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment%20smaller%20plots%20for%20larger%20text-2.png)<!-- -->

``` r
cnetplot(gse, categorySize="pvalue", showCategory = 5, foldChange=gene_list) 
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment%20smaller%20plots%20for%20larger%20text-3.png)<!-- -->

``` r
heatplot(gse, foldChange = gene_list, showCategory = 10, ) + scale_fill_gradient2(low = "red", , high = "blue", midpoint = 0)
```

![](Analysis_files/figure-gfm/Gene%20Ontology%20enrichment%20smaller%20plots%20for%20larger%20text-4.png)<!-- -->

``` r
# # Perform functional enrichment analysis on list this makes an HTML object
# gost_info <- gost(query = unique_genes)
#  
# # This information can be used to enrich these genes for further pathway analysis
# gost <- gostplot(gost_info)
```
