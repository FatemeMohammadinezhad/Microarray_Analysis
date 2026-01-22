


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
normalized_data<-rma(raw_data)
arrayQualityMetrics(expressionset = normalized_data[,c(16,17,35,36,39,41,45,46,
                                                       65,66,95,96,
                                                       99,102,106)],
                    outdir = "Results_lungcancer/QC_normalized_data",
                    force = TRUE,
                    reporttitle = "QC normalized subset (15 outlier samples)")

# In the normalized results, we can see that the problem of outlier samples now are solved # 

class(normalized_data)
processed_data<-as.data.frame(exprs(normalized_data))
dim(processed_data)
#-------------------------------------------------------------------------------

### Filter low-values probes ==> not strongly contributing in biological info
row_median<-rowMedians(as.matrix(processed_data))
row_median
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution",
     col = "blue",
     border = "white")
dev.off()
dev.list()
threshold=4.5
abline(v=threshold,col="red",lwd=2.5)
index<-row_median>threshold
filtered_data<-processed_data[index,]

### Change col names based on the phenotype data raw names
colnames(filtered_data)<-rownames(phenotype_data)
processed_data<-filtered_data
#-------------------------------------------------------------------------------
### Classifying target col names in phenotype data

class(phenotype_data$source_name_ch1)
groups<-factor(phenotype_data$source_name_ch1,
               levels = c("frozen tissue of primary tumor",
                          "frozen tissue of adjacent normal"),
               labels = c("cancer","normal"))
class(groups)
levels(groups)
