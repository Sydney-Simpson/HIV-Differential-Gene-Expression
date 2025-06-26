Data Cleaning and Preprocessing
================
2025-05-20

``` r
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

``` r
#Load required packages
library(GEOquery)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
library(limma)
```

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
library(umap)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:generics':
    ## 
    ##     explain

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(illuminaHumanv4.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: org.Hs.eg.db

    ## 

    ## 

``` r
library(AnnotationDbi)
library(tidyselect)
library(stringr)
library(sva)
```

    ## Loading required package: mgcv

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     collapse

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## This is mgcv 1.9-3. For overview type 'help("mgcv-package")'.

    ## Loading required package: genefilter

    ## Loading required package: BiocParallel

``` r
# list of GEO series IDs
geo_series <- c("GSE195434", "GSE200606")

# initialize lists for results
raw_list <- c()
sample_info_list <- c()
exprs_list <- c()

# loop through each series
for (id in geo_series) {
  # download the series
  raw_data <- getGEO(id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]

  # store raw data
  raw_list[[id]] <- raw_data
  
  # extract sample info
  sample_info_list[[id]] <- pData(raw_data)
  
  # (Optional) extract and store expression matrix. However we want to do our own normalization 
  exprs_list[[id]] <- exprs(raw_data)
}
```

    ## Found 1 file(s)

    ## GSE195434_series_matrix.txt.gz

    ## Found 1 file(s)

    ## GSE200606_series_matrix.txt.gz

``` r
#Make a new list for the annotated expression data
Annotated <- list()
  
for (set in 1:length(exprs_list)) {
  dataset <- data.frame(exprs_list[[set]])
  probe_ids <- rownames(dataset)
  dataset <- dataset %>%
    mutate(symbol = mapIds(illuminaHumanv4.db,
                       keys=probe_ids,
                       column="SYMBOL",
                       keytype="PROBEID",),
           .before = colnames(dataset)[1])
  #Remove any NAs
  dataset <- na.omit(dataset)
  #Add to the Annotated list
  Annotated[[set]] <- dataset
  # Ensure correct names
  names(Annotated)[[set]] <- names(exprs_list)[[set]]
}
```

    ## 'select()' returned 1:many mapping between keys and columns
    ## 'select()' returned 1:many mapping between keys and columns

``` r
# Select relevant columns from patient information
p_data_1 <- sample_info_list[[1]]
selected_info_1 <- p_data_1 %>%
  dplyr::select(
    geo_accession = "geo_accession",
    patient_id = "individual:ch1",
    status = "group:ch1",
    time = "time:ch1"
    ) %>%
  mutate(time = str_split_i(time, " ", 1))

# Substitute values on infection status so both data sets match
selected_info_1$status <- gsub("Post.infection","HIV-1 infected", selected_info_1$status)
selected_info_1$status <- gsub("Pre.infection","HIV negative", selected_info_1$status)


p_data_2 <- sample_info_list[[2]]
selected_info_2 <- p_data_2 %>%
   dplyr::select(
     geo_accession = "geo_accession",
     status = "source_name_ch1",
     treatment = "treatment:ch1",
     gender = "gender:ch1",
     age = "age:ch1"
  )

#save relevant sample information
write.table (selected_info_1, "GSE195434_sample_info.txt", row.names= TRUE, quote= FALSE, sep="\t")
write.table (selected_info_2, "GSE200606_sample_info.txt", row.names= TRUE, quote= FALSE, sep="\t")

#save this information as an object
sample_info <- c()
sample_info[[1]] <- selected_info_1
sample_info[[2]] <- selected_info_2

names(sample_info) <- names(sample_info_list)
```

``` r
# Make a new object to hold your processed data
p_data <- bind_rows(sample_info) %>%
  dplyr::select(geo_accession, status)

# Perform log2 transoformation on datasets and combine
d1 <- Annotated[[1]][,-1]
d2 <- Annotated[[2]][,-1]
data <- cbind(log2(d1), log2(d2))

# Perform ComBat analysis to correct for any batch effects, while preserving varaiation between conditions
batch <- c(rep( 1, times = ncol(d1)), rep(2, times = ncol(d2)))
mod <- model.matrix(~as.factor(status), data = p_data)
final <- ComBat(data, mod = mod, batch = batch, par.prior = TRUE)
```

    ## Found2batches

    ## Adjusting for3covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding parametric adjustments

    ## Adjusting the Data

``` r
# Restore gene symbol information
symbols <- data.frame("GeneSymbol" = Annotated[[1]][,1])
processed_data <- cbind(symbols, final)

# Take median value for gene symbols with multiple probes
processed_dataset <- processed_data %>%
   group_by(GeneSymbol) %>%
   summarize(across(everything(), \(x) median(x, na.rm = TRUE))) %>%
   tibble::column_to_rownames(var = "GeneSymbol")

#Save the final table
write.table(processed_dataset, "Processed_data.txt", row.names = TRUE, quote = FALSE, sep = "\t")
```

``` r
s_info_1 <- selected_info_1 %>%
  mutate(dataset = "GSE195434")

# Make sure both sets have matching columns
s_info_2 <- selected_info_2 %>%
  dplyr::select(geo_accession , status) %>%
  mutate(patient_id = c(1:nrow(selected_info_2)) ) %>%
  mutate (time = NA) %>%
  mutate(dataset = "GSE200606")

all_sample_info <- rbind(s_info_1, s_info_2)

#Save the final table
write.table(all_sample_info, "all_sample_info.txt", row.names = TRUE, quote = FALSE, sep = "\t")
```
