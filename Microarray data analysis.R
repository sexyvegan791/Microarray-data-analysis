# install packages
BiocManager::install("limma")

# upload required packages 
library(GEOquery)
library(affy)
library(dplyr)
library(limma)
library(tidyverse)
library(ggplot2)
library(data.table)


# download microarray file
getGEOSuppFiles("GSE20966")

# load count matrix file 
untar("GSE20966_RAW.tar", exdir = "GEO")
cels = list.files("GEO/", pattern = "CEL|cel")

# sometiles, it is 'CEL', you need to check it first
sapply(paste("GEO", cels, sep = "/"),gunzip)

# specifiy the file location
cels = list.files("GEO", pattern = "CEL|cel")         # specifiy the file location before running this command

# read CEL files
raw.data = ReadAffy(verbose = TRUE, filenames = cels)


# see your data with box plot before normalization
par(mar=c(9,4,3,7))
boxplot(raw.data, las=2, col=c(rep("blue",10),rep("red",10)),notch=TRUE, main="Before Normalization" )
mtext("",side = 1, line = 5, cex = 1.5)         
mtext("Expression Values", side=2, line = 2)
legend("top", c("ND", "T2D"), border="black",    # add legends
       fill = c("blue", "red"),inset = .005,
       horiz=TRUE ,cex = 0.6)



# Perform Normalization
normalized.data <- affy::rma(raw.data)


# Get normalized expression data
normalized.exprs <- data.frame(exprs(normalized.data))

# Check how your data looks like after normalization
par(mar=c(9,4,3,7))
boxplot(normalized.exprs, las=2, col=c(rep("blue",10),rep("red",10)),notch=FALSE, main="After Normalization" )
mtext("",side = 1, line = 5, cex = 1.5)
mtext("Expression Values", side=2, line = 2)
legend("topright", c("ND", "T2D"), border="black",    # add legends
       fill = c("blue","red"),inset = .005,
       horiz=TRUE ,cex = 0.6)


# Map probe IDs to Gene Symbols (Needs to download series matrix file)
gset <- getGEO("GSE20966", GSEMatrix = TRUE)


# Extract phenodata from series matrix file
phenoData <- pData(phenoData(gset[[1]]))
head(phenoData)



# Data manupulation
phenoData <- phenoData[,c(2,45)]       # choose required column
phenoData <- phenoData %>% 
  rename("diabetes_status"="disease:ch1")
#  rename("diabetes_status"="diagnosis:ch1") %>%
#  rename("age"= "age of patient:ch1") %>%
#  rename("gender"="sex:ch1")

# save file
write.csv(phenoData,"phenoData.csv")

# extract feature data
feature.data <- as.data.frame(gset[[1]]@featureData@data)

# choose required column
feature.data <- feature.data[,c(1,12)]

# combine normalized expression and feature data based on IDs (inner joining)
normalized.exprs <- normalized.exprs %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(.,feature.data, by='ID')


# export expression data to another variable 
data <- normalized.exprs

# removing ID column 
data <- dplyr::select(data,-c("ID"))

# Since there are duplicated gene symbols
# calculating means for every row grouping by the gene; and selecting rows with max mean
data <- setDT(data)[, .SD[which.max(rowMeans(.SD))], by=`Gene Symbol`]

# Removing if there are any NA's in Gene Symbols
data <- data[!(is.na(data$`Gene Symbol`)),]
data <- data[!(data$`Gene Symbol` == ''),]

# rename coumn name
data <- data %>% 
  rename("Gene_Symbol"="Gene Symbol")


# create another file to get column names from expression matrix
# if you want your columns names without .CEL format (optional)
Geo_accession <- data.frame(colnames(data))
Geo_accession <- Geo_accession[c(2:21),]
Geo_accession <- data.frame(Geo_accession)

# comnined phenoData andGeo_accession
phenoData <- cbind(phenoData,Geo_accession)

# get the expression data with gene symbols as rownames and sample id as column names

data <- data %>% 
  gather(key = 'samples', value = 'counts', -Gene_Symbol) %>% 
  inner_join(., phenoData, by = c('samples' = 'Geo_accession')) %>% 
  dplyr::select(1,3,4) %>% 
  spread(key = 'geo_accession', value = 'counts') %>% 
  column_to_rownames(var = 'Gene_Symbol')


# moving all the information stored in data back to normalized.exprs
normalized.exprs <- data


# Saving output in a csv or R file 
write.csv(normalized.exprs, file = "normalized_exprs.csv")
save(normalized.exprs, file ="normalized_exprs.RData")



# DEGs Analysis 
normalized.exprs <- get(load("normalized_exprs.RData"))
phenodata <- read.csv("phenoData.csv")



# Design and fit model
design <- model.matrix(~0+diabetes_status, phenodata)
colnames(design) <- c("ND", "T2D")
fit <- lmFit(normalized.exprs, design)


# Make contrast and fit 
cont.matrix <- makeContrasts(T2D - ND , levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)


# Get top table
tT <- topTable(fit2, adjust.method = "BH", sort.by="B", number=Inf)

# Save data 
write.csv(tT, file="DEGs.csv", row.names = T)
save(tT, file="DEGs.RDA")



# Screen over-expressed and under-expressed genes based on adj.Pval and logFC
DEGs <- read.csv("DEGs.csv")

normalized.exprs <- read.csv("normalized_exprs.csv")
up <- subset(DEGs, logFC > 0.5 & P.Value < 0.05)
down <- subset(DEGs, logFC < -0.5 & P.Value < 0.05)
