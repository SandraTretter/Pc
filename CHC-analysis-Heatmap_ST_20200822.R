# Sandra Tretter
# 19.08.2020


#### Pogonomyrmex californicus, CHC analysis, Heatmap

#---------------------------------------------------------------------

# content:  0. required packages/libraries

#           1. read in and transform the data

#           2. Heatmap

#---------------------------------------------------------------------

#### 0. required packages/libraries ####

# load required packages
library(afex)
library(MASS)
library(car)
library(tidyr)
library(dplyr)
library(broom)
library(emmeans)
#devtools::install_github("tomwenseleers/export")
library(export)
library(pheatmap)
library(doBy)
library(multcomp)
library(ggplot2)
library(ggthemes)
library(scales)
library(gtools)
library(matrixStats)
library(ggpubr)


#---------------------------------------------------------------------

#### 1. read in and transform the data ####

# set your working directory
setwd("C:/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten")

# read in your data
batch = read.csv("/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten/BatchTable_CSV_ST_20200814.csv",check.names=F,stringsAsFactors=FALSE, sep=";",dec=",") # raw peak areas

# create a subset of the sample names to add this later to your standardizied values 
samplenames <- subset(batch, select = "Sample")  

# create a subset of only your values for standardization because you cannot standardize your data with non-numeric values (Sample names)
batchvalues <- subset(batch, select = -c(1)) 

# get rid of NAs
batchvalues[is.na(batchvalues)] = 0


### 2 different standardizations:

## 1. standardization with internal C12-standard to compensate deviations caused by injection volume

for (i in 1:nrow(batchvalues)) {
  
  batchvalues[i,]	<-  batchvalues[i,] / batchvalues[i,1]
  
}

## 2. relative responses to make peaks between different samples comparable 

# if standard is the first coloumn in your dataset
for (i in 1:nrow(batchvalues)) {
  
  batchvalues[i,]	<- batchvalues[i,] / sum(batchvalues[i,])
  
}

## if standard is not the first column
#for (i in 1:nrow(batch)) {
#  batch[,2:ncol(batch)]	<- batch[,2:ncol(batch)] / sum(batch[,2:ncol(batch)])
#}

# get rid of the internal standrad column
batchvalues = subset(batchvalues, select = -c(1))  

# put your sample names and standardized values in one dataframe together
batch <- cbind(samplenames, batchvalues)  


batch[,2:ncol(batch)] = as.numeric(as.matrix(batch[,2:ncol(batch)]))

rownames(batch) = make.unique(as.character(batch[,1]))


full_batch = batch # dataframe with factor Sample + abundances
batch_factors = batch[,1,drop=FALSE] # dataframe with factor Sample only
batch = batch[,-1] # dataframe with abundances only
t_batch = t(batch) # same but transposed

batch_log10_t = t(log10(batch+1E-6)) # log10 transformed relative peak areas (transposed)

# we also calculate row z scores on log2 transformed relative peak areas for heatmap :
batch_log10_t_centered = sweep(batch_log10_t,1,rowMeans(batch_log10_t),FUN="-") 
batch_log10_t_zscores = sweep(batch_log10_t_centered,1,rowSds(batch_log10_t_centered),FUN="/") 


#---------------------------------------------------------------------

#### 2. Heatmap ####

# Heatmap of row z scores calculated on log2 transformed relative peak areas ####

 # we clip heatmap scale between -2 and 2
batch_log10_t_zscores[batch_log10_t_zscores>2] = 2  # everything above 2 is 2
batch_log10_t_zscores[batch_log10_t_zscores<(-2)] = -2  # everything under 2 is -2
ph=pheatmap(batch_log10_t_zscores, 
            cluster_rows = FALSE,  # no cluster for the rows, if TRUE than also cutree_rows for how many rows 
            cutree_cols=3,  # in how many "columns" should the heatmap be divided
            clustering_method="average", 
            clustering_distance_cols="correlation",
            #cellheight = 10,
            #cellwidth = 12
            )
#comporder_heatmap=rownames(batch_log10_t_zscores)[ph$tree_row$order]
graph2png(file="Fig1_heatmap.png",width=12,height=21,dpi=600) # in Fig. 1 in MS Q & VQ groups were manually switched around
graph2ppt(file="Fig1_heatmap.pptx",width=12,height=21) 
