# Microarray-data-analysis
# The download all the packages required for the analysis 
# R code to download Bioconductor packages 
BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("oligo") # This is optional

# install R packages for data manipulation
install.packages("dplyr")
install.packages("tidyverse")
install.packages("data.table")

# Now these packages needs to import into R workplace
library(GEOquery)
library(affy)
library(limma)
library(oligo)
library(dplyr)
library(tidyverse)
library(data.table)

# stats with the analysis followed by several steps
