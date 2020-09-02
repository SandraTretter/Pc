# Sandra Tretter
# 26.08.2020


#### Pogonomyrmex californicus, CHC Analyse

#---------------------------------------------------------------------

# content:  0. required packages/libraries

#           1. read in and transform the data

#           2. CHC classes

#           3. Barplots

#---------------------------------------------------------------------

#### 0. required packages/libraries ####

library(permute)
library(lattice)
library(vegan)
library(MASS)
library(gtools)
library(shape) 
library(ggplot2)
library(viridis)
library(export)

#---------------------------------------------------------------------

#### 1. read in and transform the data ####

setwd("C:/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten")

dataset <- read.csv("/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten/BatchTable_CSV_ST_20200814.csv" ,header=T,dec=",",sep=";",check.names=FALSE,row.names = 1)


# Get rid of NAs

dataset[is.na(dataset)] = 0


### 2 different standardizations:

## 1. standardization with internal C12-standard to compensate deviations caused by injection volume

for (i in 1:nrow(dataset)) {
  
  dataset[i,]	<- dataset[i,]*22.5 / dataset[i,1]    # internal standard is here the first column 
  
}

## 2. relative responses to make peaks between different samples comparable 

for (i in 1:nrow(dataset)) {
  
  dataset[i,]	<- dataset[i,] / sum(dataset[i,])
  
}


# Get rid of the internal standard column

dataset = subset(dataset, select = -c(1))

#-----------------------------------------------------------------------

#### 2. CHC classes ####

nm <- colnames(dataset)

alkanes <- nm %in% grep("n-", nm, value = TRUE)
alkanes <- subset(dataset, select = alkanes)
alkanes$sample<-row.names(alkanes)  # add column with sample names to data frame for plotting them on the x axis of the barplots 
alkanes$caste<-gsub("Pc_(.*?)_(.*?)_.*","\\1",alkanes$sample,perl=T)  # add column with caste to data frame for plotting them on x axis of boxplots
alkanes$lineage<-gsub("Pc_(.*?)_(.*?)_.*","\\2",alkanes$sample,perl=T)  # add lineage to data frame 
df<-tidyr::gather(alkanes,key="chc",value="val",-sample,-caste,-lineage)
dfAlkanes <- df
dfAlkanes$group <- "alkanes"

alkenes <- nm %in% grep("alkene", nm, value = TRUE)
alkenes <- subset(dataset, select = alkenes)
alkenes$sample<-row.names(alkenes)
alkenes$caste<-gsub("Pc_(.*?)_(.*?)_.*","\\1",alkenes$sample,perl=T)
alkenes$lineage<-gsub("Pc_(.*?)_(.*?)_.*","\\2",alkenes$sample,perl=T)
df<-tidyr::gather(alkenes,key="chc",value="val",-sample,-caste,-lineage)
dfAlkenes <- df
dfAlkenes$group <- "alkenes"

me <- nm %in% grep("-Me", nm, value = TRUE)
me <- subset(dataset, select = me)
me$sample<-row.names(me)
me$caste<-gsub("Pc_(.*?)_(.*?)_.*","\\1",me$sample,perl=T)
me$lineage<-gsub("Pc_(.*?)_(.*?)_.*","\\2",me$sample,perl=T)
df<-tidyr::gather(me,key="chc",value="val",-sample,-caste,-lineage)
dfMe <- df
dfMe$group <- "monomethyl-alkanes"

di <- nm %in% grep("-Di", nm, value = TRUE)
di <- subset(dataset, select = di)
di$sample<-row.names(di)
di$caste<-gsub("Pc_(.*?)_(.*?)_.*","\\1",di$sample,perl=T)
di$lineage<-gsub("Pc_(.*?)_(.*?)_.*","\\2",di$sample,perl=T)
df<-tidyr::gather(di,key="chc",value="val",-sample,-caste,-lineage)
dfDi <- df
dfDi$group <- "dimethyl-alkanes"

tri <- nm %in% grep("-Tri", nm, value = TRUE)
tri <- subset(dataset, select = tri)
tri$sample<-row.names(tri)
tri$caste<-gsub("Pc_(.*?)_(.*?)_.*","\\1",tri$sample,perl=T)
tri$lineage<-gsub("Pc_(.*?)_(.*?)_.*","\\2",tri$sample,perl=T)
df<-tidyr::gather(tri,key="chc",value="val",-sample,-caste,-lineage)
dfTri <- df
dfTri$group <- "trimethyl-alkanes"

tetra <- nm %in% grep("-Tetra", nm, value = TRUE)
tetra <- subset(dataset, select = tetra)
tetra$sample<-row.names(tetra)
tetra$caste<-gsub("Pc_(.*?)_(.*?)_.*","\\1",tetra$sample,perl=T)
tetra$lineage<-gsub("Pc_(.*?)_(.*?)_.*","\\2",tetra$sample,perl=T)
df<-tidyr::gather(tetra,key="chc",value="val",-sample,-caste,-lineage)
dfTetra <- df
dfTetra$group <- "tetramethyl-alkanes"

diene <- nm %in% grep("diene", nm, value = TRUE)
diene <- subset(dataset, select = diene)
diene$sample<-row.names(diene)
diene$caste<-gsub("Pc_(.*?)_(.*?)_.*","\\1",diene$sample,perl=T)
diene$lineage<-gsub("Pc_(.*?)_(.*?)_.*","\\2",diene$sample,perl=T)
df<-tidyr::gather(diene,key="chc",value="val",-sample,-caste,-lineage)
dfDiene <- df
dfDiene$group <- "dienes"


# combine classes
dfall <- rbind(dfAlkanes, dfAlkenes, dfMe, dfDi, dfTri, dfTetra, dfDiene)

dfall$group <- factor(dfall$group, 
                      levels = c("alkanes","alkenes","dienes","monomethyl-alkanes","dimethyl-alkanes","trimethyl-alkanes","tetramethyl-alkanes")) # sorting bars

#-----------------------------------------------------------------------

#### 3. Barplots ####


p <- ggplot(dfall, aes(fill=group, y=val, x=caste)) +
  geom_bar(position = "fill", stat="identity") + # position = fill -> norming to 1
  theme_classic(base_size = 22) +  # white background
  theme(axis.title.x = element_blank(), axis.title.y = element_text(color = "black", size=14, face="bold")) +
  theme(legend.title =element_blank(), legend.text = element_text(color = "black", size = 16)) +
  scale_x_discrete(guide= guide_axis(n.dodge = 2)) +
  NULL +
  scale_fill_manual(values = c("goldenrod1", "darkorange2", "chartreuse3", "slategray1", "skyblue2", "steelblue3", "steelblue4")) +
  facet_grid(.~lineage, scales = "free")

p

graph2ppt(p, file="Barplots-abs_ST_20200829",width=15,height=11)


# calculate relative values to plot in bars

relValues <- dfall %>% group_by(caste, lineage, group) %>% summarise(val = sum(val))   # sum of all values per group
malehaploall <- dfall %>% group_by(caste, lineage) %>% summarise(val = sum(val))    # sum of all values --> sum of val per group / sum of all val per group -> for every caste/lineage combi

together <- dfall %>% group_by(caste, lineage, group) %>% summarise(val = sum(val)/malehaploall[1,2])
