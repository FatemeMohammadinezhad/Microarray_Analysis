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
ls("package:hgu133plus2.db")
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
  arrange(desc(probes_per_gene))

duplicate_summary

#identify genes associated with multiple probes
duplicate_genes<-duplicate_summary%>%
  filter(probes_per_gene>1)

sum(duplicate_genes$probes_per_gene)
################################################################################

### Merge annotation mapping with expression data ###

################################################################################
# Verify if Probe IDs in mapping correspond to expression data
all(gene_map_df$PROBEID==row.names(processed_data))

# merge annotation(symbol) with expression matrix
processed_data_df<-processed_data%>%
  as.data.frame()%>%
  tibble::rownames_to_column("PROBID")%>%
  dplyr::mutate(SYMBOL=gene_symbols[PROBID])%>%
  dplyr::relocate(SYMBOL,.after=PROBID)

#Remove prob without valid gene symbols annotation
processed_data_df<-processed_data_df%>%
  dplyr::filter(!is.na(SYMBOL))

#Select only numeric expression column
expr_only<-processed_data_df%>%
  dplyr::select(-PROBID,-SYMBOL)
#-------------------------------------------------------------------------------

### limma::avereps() computes the average for probes that represent or multiple same genes
# Example of avereps work
x<-matrix(rnorm(8*3),8,3)
colnames(x)<-c("S1","S2","S3")
rownames(x)<-c("b","a","a","b","c","c","a","b")
head(x)
avereps(x)  #collapse duplicated renames by averaging
#-------------------------------------------------------------------------------

averaged_data<-limma:: avereps(expr_only,ID=processed_data_df$SYMBOL)
dim(averaged_data)

data<-as.data.frame(averaged_data)
data<-data.matrix(data)
str(data) #check structure data
is.numeric(data)

#-------------------------------------------------------------------------------
#### Differential Gene Expression Analysis ####
#-------------------------------------------------------------------------------
# Define sample groups based on phenotype data
# Adjust group labels according to data set annotation

groups<-factor(phenotype_data$source_name_ch1,
               levels = c("gastric mucosa","gastric adenocarcinoma"),
               labels = c("normal","cancer"))

class(groups)
levels(groups)
#-------------------------------------------------------------------------------
# Create design matrix for linear modeling 
#-------------------------------------------------------------------------------
# Using no intercept(~0+groups) allows each group to have its own coefficient
design<-model.matrix(~0+groups)
colnames(design)<-levels(groups)

# Fit linear model to expression data
fit_1 <- lmFit(data, design)


# Define contrast to compare cancer vs normal samples
contrast_matrix <- makeContrasts(cancer_vs_normal = cancer - normal,
                                 levels = design)

# Apply contrasts and compute moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# -------------------------------------------------------------
# Extract list of differentially expressed genes (DEGs)
# -------------------------------------------------------------
# -------------------------------------------------------------
# Extract list of differentially expressed genes (DEGs)
# -------------------------------------------------------------
deg_results <- topTable(fit_2,
                        coef = "cancer_vs_normal",  # Specify contrast of interest
                        number = Inf,               # Return all genes
                        adjust.method = "BH")       # Benjamini-Hochberg correction

# -------------------------------------------------------------
# Classify DEGs into Upregulated, Downregulated, or Not Significant
# -------------------------------------------------------------
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))
# Subset genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")
# Combine both sets of DEGs
deg_updown <- rbind(upregulated, downregulated)

write.csv(deg_results, file = "Results/DEGs_Results.csv")
write.csv(upregulated, file = "Results/Upregulated_DEGs.csv")
write.csv(downregulated, file = "Results/Downregulated_DEGs.csv")
write.csv(deg_updown, file = "Results/Updown_DEGs.csv")
# -------------------------------------------------------------
#### Data Visualization ####
# -------------------------------------------------------------

# -------------------------------------------------------------
# Volcano Plot: visualizes DEGs by logFC and adjusted p-values
# -------------------------------------------------------------
# Note: x-axis = log2 fold change, y-axis = -log10 adjusted p-value

# Save volcano plot as PNG

dir.create("Result_Plots", showWarnings = FALSE)

png("Result_Plots/volcano_plot.png", width = 2000, height = 1500, res = 300)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(P-value)",
       color = "Regulation")

dev.off()

# -------------------------------------------------------------
# Heatmap of Top Differentially Expressed Genes
# -------------------------------------------------------------

# Select top genes with smallest adjusted p-values
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 10)

# Subset averaged expression matrix for selected genes
heatmap_data <- data[top_genes, ]

# Generate unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

# Assign formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names



dev.off()

dev.cur()
dev.off()

png("heatmap_top10_DEGs.png", width = 2000, height = 1500, res = 300)

pheatmap(
  heatmap_data,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 10 Differentially Expressed Genes"
)

dev.off()
file.info("heatmap_top10_DEGs.png")$size



