# Clean start
rm(list = ls())

# Ensure BiocManager exists
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}

# Set Bioconductor version for R 4.5.x
BiocManager::install(version = "3.22", ask = FALSE)

# IMPORTANT (Windows): prefer binaries so you DON'T need Rtools
options(pkgType = "binary")

# Make sure Bioconductor repos are correctly set
options(repos = BiocManager::repositories())

# Install dependencies + the package (force reinstall)
BiocManager::install(c(
  "arrayQualityMetrics",
  "GEOquery",
  "affy",
  "limma",
  "AnnotationDbi",
  "hgu133plus2.db"
), ask = FALSE, update = TRUE, dependencies = TRUE, force = TRUE)

install.packages(c("dplyr","tibble","ggplot2","pheatmap"), dependencies = TRUE)
################################################################################
library(arrayQualityMetrics)
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(limma) # perform linear modeling & differential expression
library(AnnotationDbi)
library(hgu133plus2.db) # annotation matrix for affymatrix
library(dplyr) #simplify data manipulation tasks # belongs to tidyverse
library(tibble)
library(ggplot2) #belongs to tidyverse
library(pheatmap)

################################################################################
### Download series Matrix files ####
################################################################################

gse_data<-getGEO("GSE79973",GSEMatrix = TRUE)
#GSEMatrix = TRUE, tells the function to return a series matrix object which contains 
#process expression data along with phenotype & feature data.

# Extract expression data
expression_data<-exprs(gse_data$GSE79973_series_matrix.txt.gz)

# Extract feature data
feature_data<- fData(gse_data$GSE79973_series_matrix.txt.gz)

# Extract phenotype data
phenotype_data<-pData(gse_data$GSE79973_series_matrix.txt.gz)
sum(is.na(phenotype_data$source_name_ch1)) 

################################################################################
### Download raw data(CEL files)
################################################################################

getGEOSuppFiles("GSE79973",baseDir = "R_Projects_Raw_Files",makeDirectory = TRUE)

### Untar CEL files
untar("R_Projects_Raw_Files/GSE79973/GSE79973_RAW.tar",exdir ="R_Projects_Raw_Files/CEL_files" )

### read CEL files
raw_data<-ReadAffy(celfile.path = "R_Projects_Raw_Files/CEL_files")
raw_data
### NOTE:If the data set is coming from "ilumina" or "agilant",we may need 
### additional bioconductor packages to read and pre-process the data
### EXAMPLE: lumi package for ilumina

################################################################################
### QC before pre-processing data
# box plots
# MA plots
# density plots

### This array quality metrics packages is specified for micro array data
################################################################################

arrayQualityMetrics(expressionset =raw_data,
                    outdir = "Results/QC_raw_data",
                    force = TRUE,
                    do.logtransform = TRUE )

################################################################################
### RMA Normalization
################################################################################
### 
#1️⃣ Background correction==>Removes optical noise and nonspecific binding.

#2️⃣ Quantile normalization==>Forces all samples to have the same overall distribution, making them comparable.

#3️⃣ Probe stigmatization==>Each gene has many probes → rma() combines them into one expression value per gene using a robust statistical mod

normalized_data<-rma(raw_data)

### QC after normalization
arrayQualityMetrics(expressionset =normalized_data[,c(1,8,9,11,20)],#bcz these samples were outlier
                    outdir = "Results/QC_normalized_data",
                    force = TRUE,
                   )

### Extract normalized expression values
processed_data<-as.data.frame(exprs(normalized_data))
dim(processed_data)

### Now we want to remove probes that show little variation across the samples 
### so we filter low values probes, bcz probes with constant expression  across 
### all samples usually do not contribute in biological information or it may 
### introduce noise ==> do median of prob sets for each gene .
row_median<-rowMedians(as.matrix(processed_data))
row_median

### plot distribution of median intensities of prob sets
hist(row_median,breaks=100,freq=FALSE, main="Median Intensity Distribution",
     col="green",border="white")

threshold<-3.5
abline(v=threshold,col="red",lwd = 2)
indx<-row_median>threshold
indx

filtered_data<-processed_data[indx, ]

### change the long column names based on the phenotype_data names
colnames(filtered_data)<-row.names(phenotype_data)
processed_data<-filtered_data

################################################################################
### Phenotype data
################################################################################

class(phenotype_data$source_name_ch1)
groups<-factor(phenotype_data$source_name_ch1,
               levels = c("gastric mucosa","gastric adenocarcinoma"),
                          label=c("normal","cancer"))
class(groups)
levels(groups)
################################################################################
####Prob_to_gene Mapping using AnnotationDBI####

annotation(raw_data)
raw_data

ls("package:hgu133plus2.db")

columns(hgu133plus2.db)
keytypes(hgu133plus2.db)
probe_ids<-rownames(processed_data)
gene_symbols<-mapIds(hgu133plus2.db,
                     keys=probe_ids,
                     keytype="PROBEID",
                     column="SYMBOL",
                     multiVals = "first")

symbols<-AnnotationDbi::select(hgu133plus2.db,
                               keys=probe_ids,
                               keytype="PROBEID",
                               columns =c("SYMBOL","ENTREZID","GENENAME"))

gene_map_df<-gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL=2)

### Handle multiple probes mapping to a single gene

### Several strategies:
#1.Retain probes with highest expression or variance 
#2.Average or summarize probes signals
#3.Remove duplicate probes to maintain single probe for each gene

#summarize probes signals
duplicate_summary<-gene_map_df%>%
  group_by(SYMBOL)%>%
  summarise(probes_per_gene=n())%>%
  arrange(desc(c))

duplicate_summary

#identify genes associated with multiple probes
duplicate_genes<-duplicate_summary%>%
  filter(probes_per_gene>1)
sum(duplicate_genes$probes_per_gene)
################################################################################

### Merge annotation mapping with expression data ###

################################################################################

