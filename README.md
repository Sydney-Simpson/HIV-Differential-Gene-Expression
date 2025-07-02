# HIV Differential Gene Expression
# Table of Contents
- [Background](#Background)
- [Preprocessing and Cleaning](#Preprocessing)
- [Differential Expression](#Differential-Expression)
- [Principal Component Analysis](#Principal-Component-Analysis)
- [Gene Enrichment Analysis Using Gene Ontology](#Gene-Enrichment-Analysis-Using-Gene-Ontology)

# Background
Use of microarray data for gene expression analysis is widely used to study changes in biological processes across different disease states. Whole microarrays consists of tens of thousands of probes corresponding to different genes.
The steps in analyzing microarray data include quality control, cleaning, normalization, differential expression analysis and potential downstream applications. This project combines datasets from two separate studies and thus includes a step to correct for batch effects.

This project aims to perform a microarray gene expression analysis on samples from HIV positive and negative persons generated from two independent studies.

## Objective
The objective of this project is to identify genes which are differentially expressed in HIV positive and HIV negative persons.

## Dataset Descriptions
Two datasets from the GeneExpression Omnibus (GEO) are used for this analysis. Both are sets of microarrays containing gene expression data from samples of blood from HIV positive and negative persons.
Both datasets utilized the Illumina HumanHT-12 V4.0 expression beadchip as a platform. The combined dataset consisted of 26 HIV negative and 45 HIV positive samples. 
HIV positive samples were further broken down into 30 HIV-1 positive, 7 HIV-2 positive, and 8 HIV-1 and HIV-2 dual positive persons.

[GSE195434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195434) <br/>
 This study included microarray expression data from 26 participants, with samples taken periodically over two years. 18 subjects were HIV negative, but later became HIV positive during the course of the study. For this analysis the last sample while individuals were 
 still HIV negative, and the first sample following a positive test were used. The first results from individuals who entered the study with an HIV positive status were employed for the other 8 participants. <br/>
 <br/>
[GSE200606](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200606) <br/>
Enrolled 31 participants in total, including HIV negative individuals as well as persons who tested positive for HIV-1, HIV-2, or both HIV-1 and HIV-2. There 8 participants per group for all status conditions except for the group of people infected with HIV-2, 
which only included 7 participants. <br/>

Note: That in this file either R or a SQL script (using sqldf) can be used for sample selection.
# Preprocessing
All data was obtained from GEO utilizing the GEOquery package. Expression arrays and descriptive information for all samples were extracted from both datasets. There were initially 47,231 probes or datapoints for each sample. Control probes were removed, and all probe IDs were matched to the corresponding gene symbol. Median values were used for genes with multiple corresponding probes, which resulted in 21,186 genes. <br/>
Both datasets had already undergone quantile normalization, they were further log2 transformed then combined using ComBat to minimize batch effects. Later correlation plots showed little to no grouping by original dataset. 

![create heatmaps and narrow down genes list-3](https://github.com/user-attachments/assets/7ce05b7f-7f24-4a79-98c9-79bbadaba02b)

# Differential Expression
Mean fold change values were calculated for each HIV positive samples group. A linear regression was used to identify genes with statistically significant differential expression between each HIV positive group and the negative group. A Benjamini-Hochberg procedure was used to account for the false discovery rate and p values were adjusted accordingly. 120 genes showed statistically signifcant difference between samples from HIV-1 infected and HIV negative persons. Additionally 1 gene in the HIV-1 and HIV-2 dual infected population showed a statistically significant difference. Both HIV-2 and dual infected populations had a low samples sizes, this impacted these population's adjusted p value scores.

![Analyze different groups using limma-1](https://github.com/user-attachments/assets/edaada6e-3361-4de4-bf7f-9773ca189344)

The volcano plots below show the genes with the highest differences in expresion between HIV negative and HIV-1 positive persons.
|  Genes with highest fold change  |  Most statistically significant genes |
|---------|--------|
| ![Identify changes in gene expression then match with p values for each infection status-8](https://github.com/user-attachments/assets/b612be05-860e-4dff-b71f-73877c58b2ef) | ![Identify changes in gene expression then match with p values for each infection status-7](https://github.com/user-attachments/assets/d292adc6-dbd6-4fd1-9eda-e2a5c1276db0)

TMEM119 was the only gene to show statistically different significance in the dual infected population. 

![Identify changes in gene expression then match with p values for each infection status-9](https://github.com/user-attachments/assets/31fc0021-28ad-4692-b79a-467cf9d2dac3)

From here genes that were selected which either showed statistical signifcance between HIV infected and uninfected populations, or had at least a 1.5x (or log2 0.6) fold change. The resulting 225 genes were assessed further. Outliers were assessed by collecting outliers for each gene with regard to infection status. Samples with 40 or more outliers out of the 225 total genes were removed from the analysis. <br/>

Data was transformed by subtracting the average value for each individual gene. All heatplots show log2 change as compared to the average gene expression for that gene in the dataset. 

![create heatmaps and narrow down genes list-4](https://github.com/user-attachments/assets/65038404-40d6-48dc-b86e-5564485bacd7)

| Top 10 Genes with the highest average fold change  | 10 Most statistically significant genes | 
|----|----|
| ![create heat maps-4](https://github.com/user-attachments/assets/7280bd98-03ee-4957-8944-fdf747a4cc5b)| ![create heat maps-3](https://github.com/user-attachments/assets/59b3b2a4-5675-4f2f-9d91-245204b47688) |

Samples primarily clustered by HIV positive and HIV negative infection status with little to no difference observed between HIV-1, HIV-2 or Dual infected samples. Notably IF27 showed significant upregulation in HIV infected individuals, this is consistent with previous studies. Other notably genes included interferons and interferon stimulated genes (ISGs) such as IFI44 it's ligand (IFI44). Finally, 20 genes which were statistically significant and showed  an average of at least 1 fold change between HIV positive and negative samples were selected for PCA analysis.

# Principal Component Analysis
| Groups | Genes with greatest representation |
|---|----|
| ![PCA for genes with p value less than or equal to 0 05 and a 1 or greater fold change-1](https://github.com/user-attachments/assets/71d611c0-067e-4333-912b-f56dca87d2cd) | ![PCA for genes with p value less than or equal to 0 05 and a 1 or greater fold change-2](https://github.com/user-attachments/assets/4723f993-f9cb-483b-9874-226d294a8fa6) |

The elbow method determined that 3 components was the ideal number to keep, these also explain 90% of variation observed in dataset as shown in the screeplot below.

![PCA for genes with p value less than or equal to 0 05 and a 1 or greater fold change-3](https://github.com/user-attachments/assets/b7d6a99b-1394-4648-ae0e-c1a7cd2f170c)

The top 5% of loadings were chosen for each PC1. Interestingly PC1, which by itself explains over 70% of the variation, had multiple genes contributing similarly. PC2 and PC3 however had components with positive and negative contributions. 

| | Top 5% of Variables for Each Component  |
|---|-----|
| PC1 | ![PCA for genes with p value less than or equal to 0 05 and a 1 or greater fold change-6](https://github.com/user-attachments/assets/2ce1c668-97b6-41cc-9409-67cfadf23ae6) |
| PC2 | ![PCA for genes with p value less than or equal to 0 05 and a 1 or greater fold change-7](https://github.com/user-attachments/assets/c6498726-a8bd-4236-838a-a96a7cb9e770) | 
| PC3 | ![PCA for genes with p value less than or equal to 0 05 and a 1 or greater fold change-8](https://github.com/user-attachments/assets/80df1a3e-0037-425b-9026-ddd3a3d7cef1) |

The following genes were among the top 5% for at least one of the top three PCs: EPSTI1, IFI44, IFI44L, LY6E, XAF1, CD52, EOMES, GZMH, IFI27, HERC5, IFIT1, RSAD2, NKG7

![select top genes from PCA and create heatmap-1](https://github.com/user-attachments/assets/c39b1a88-3ec2-4430-8f7c-d171d5ad1177)

# Gene Enrichment Analysis Using Gene Ontology
 Gene Ontology was used to enrich the dataset of the original 225 selected genes. Due to small sample sizes for HIV-2 and dual infected persons populations, only HIV-1 positive and HIV negative samples were included for this analysis.
 First, as expected the top 10 categories were heavily related to immune system terms. The top genes were heavily interwoven within these related terms. Interestingly, negative resulation of catabolic processes was also one of the top terms represented within this gene set. Decreasing these negative regulators of processes which break down substances within cells means that there is an increase in these activities. Early on in infection HIV is known to require one such process, autophagy for viral replication. For example BCL2L1 is from the BCL2 family of proteins, which regulate apoptosis (programmed cell death) and autophagy. BCL2L1 can be alternatively spliced into BCL2, which is known be regulated by viral proteins during infection. I actually wrote a [review article](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.682624/full) on how HIV promotes BCL2 ubinquination while I was a postdoc!

![Gene Ontology enrichment smaller plots for larger text-1](https://github.com/user-attachments/assets/cbd06f8b-1216-4019-b944-999ef95a4812)


![Gene Ontology enrichment smaller plots for larger text-2](https://github.com/user-attachments/assets/a2a05b36-48d1-47f0-9945-dfcba2e59a3f)

![Gene Ontology enrichment smaller plots for larger text-3](https://github.com/user-attachments/assets/e1cf8ad3-0909-4d78-bab5-f32024d88be4)

![Gene Ontology enrichment smaller plots for larger text-4](https://github.com/user-attachments/assets/37cfa3be-7fdc-4835-bdbf-53f1daaf924c)






