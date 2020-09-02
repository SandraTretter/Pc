# Sandra Tretter
# 19.08.2020


#### Pogonomyrmex californicus, CHC analysis, PCA, colonies

#---------------------------------------------------------------------

# content:  0. required packages/libraries

#           1. read in and transform the data

#           2. PCA

#---------------------------------------------------------------------

#### 0. required packages/libraries ####

library(ggplot2)
#install.packages("viridisLite")
library(viridis)
library(permute)
library(lattice)
library(vegan)
library(MASS)
library(gtools)
library(ggbiplot)
library(zCompositions)
library(compositions)
library(export)


#---------------------------------------------------------------------

#### 1. read in and transform the data ####

# set your working directory
setwd("C:/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten")

# read in your dataset
dataset <- read.csv("/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten/BatchTable-Colonies_CSV_ST_20200814.csv" ,header=T,dec=",",sep=";",check.names=FALSE,row.names = 1)


# Get rid of NAs
dataset[is.na(dataset)] = 0


# create a vector for grouping your samples in the different colonies (which we will use for the different colours)
vec <- vector(length=length(row.names(dataset))) 
for (i in 1:length(vec)) vec[i] = substr(row.names(dataset)[i],11,13) # Depends on sample names -> e.g. Pc_haplo_Q_1_01 has 15 characters -> 4,10 = Q_haplo

# create a vector four grouping your samples in different classes, like caste and population
class <- vector(length=length(row.names(dataset))) 
for (i in 1:length(class)) class[i] = substr(row.names(dataset)[i],4,10)

# take a look at your vectors
unique(vec)
unique (class)


### 2 different standardizations:

## 1. standardization with internal C12-standard to compensate deviations caused by injection volume

for (i in 1:nrow(dataset)) {
  
  dataset[i,]	<- dataset[i,] / dataset[i,1]
  
}

## 2. relative responses to make peaks between different samples comparable 

for (i in 1:nrow(dataset)) {
  
  dataset[i,]	<- dataset[i,] / sum(dataset[i,])
  
}

# Get rid of the internal standard column
dataset = subset(dataset, select = -c(1))

#-----------------------------------------------------------------------

#### 2. PCA  ####

# centered log-ratios transformation
clr.output <- clr(dataset)
clr.output  # take a look at your transformed data

dataset.pca <- prcomp(clr.output, scale.=T, center=T) # principal component analysis -> results saved as dataset.pca
summary(dataset.pca)

unique(colony)

all.clr <- ggbiplot(dataset.pca, obs.scale=1, var.scale=1, choices=c(1,2), # choices -> which PCAs to plot
                var.axes=F, # remove arrows
                groups = class) + # groups -> which groups should be plotted (vector you created before)
  
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  
  theme_classic() + # theme_classic() for white background
  
  theme(legend.title = element_blank()) + # removes legend title
  
  scale_color_manual(name=vec, values = c("darkred", "gold", "darksalmon", "deeppink", "azure4", "forestgreen", "orange", "purple4", "black", "black", "black", "black", "black", "black")) +
  
  scale_shape_manual(name=class, values=c(15, 0, 17, 2, 19, 1)) + # pch symbols for your groups 
  
  geom_point(aes(colour=vec, shape = class, size = 3)) # creates the symbols for the groups in your defined colours
  #geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2) # if you like to add sample names

all.clr  # take a look at your PCA

#export graph
graph2png(file="PCA_colonies_ST_20200824.png",width=15,height=11,dpi=600)
graph2ppt(file="PCA_colonies_ST_20200824.ppt",width = 13, height = 10)

#-----------------------------------------------------------------------