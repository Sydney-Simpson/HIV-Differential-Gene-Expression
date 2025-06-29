# HIV Differential Gene Expression
# Table of Contents
- [Background](#Background)
- [Preprocessing and Cleaning](#Preprocessing)
- [Differential Expression](#DifferentialExpression)

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
which only included 7 participants.

# Preprocessing
All data was obtained from GEO utilizing the GEOquery package. Expression arrays and descriptive information for all samples were extracted from both datasets. There were initially 47,231 probes or datapoints for each sample. Control probes were removed, and all probe IDs were matched to the corresponding gene symbol. Median values were used for genes with multiple corresponding probes, which resulted in 21,186 genes. <br/>
Both datasets had already undergone quantile normalization, they were further log2 transformed then combined using ComBat to minimize batch effects.

# Differential Expression
Mean fold change values were calculated for each HIV positive samples group. A linear regression was used to identify genes with statistically significant differential expression between each HIV positive group and the negative group. A Benjamini-Hochberg procedure was used to account for the false discovery rate and p values were adjusted accordingly.


![Analyze different groups using limma-1](https://github.com/user-attachments/assets/edaada6e-3361-4de4-bf7f-9773ca189344)


![create heatmaps and narrow down genes list-2](https://github.com/user-attachments/assets/ebe81f34-448d-4f2e-90a5-1b5d9f29b6cd)
