


library(arrayQualityMetrics)
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(limma) 
library(AnnotationDbi)
library(hgu133plus2.db) 
library(dplyr) 
library(tibble)
library(ggplot2) 
library(pheatmap)

getwd()
#-------------------------------------------------------------------------------
### Data Retrieval ###

gse_data<-getGEO("GSE19804",GSEMatrix = TRUE)
expression_data<-exprs(gse_data$GSE19804_series_matrix.txt.gz)
dim(expression_data)
feature_data<-fData(gse_data$GSE19804_series_matrix.txt.gz)
phenotype_data<-pData(gse_data$GSE19804_series_matrix.txt.gz)
sum(is.na(phenotype_data$source_name_ch1))

getGEOSuppFiles("GSE19804",baseDir = "R_Projects_Raw_Files",
                makeDirectory = TRUE)

untar("R_Projects_Raw_Files/GSE19804/GSE19804_Raw.tar",
      exdir = "R_Projects_Raw_Files/CEL_Files_lungcancer")

raw_data<-ReadAffy(celfile.path = "R_Projects_Raw_Files/CEL_Files_lungcancer")
#-------------------------------------------------------------------------------

### Quality Control Before pre-processing the data ###

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results_lungcancer/QC_raw_data",
                    force = TRUE,
                    do.logtransform = TRUE)


