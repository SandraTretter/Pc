# Sandra Tretter
# 26.08.2020


#### Pogonomyrmex californicus, CHC Analyse

#---------------------------------------------------------------------

# content:  0. required packages/libraries

#           1. read in and transform the data

#           2. CHC classes

#           3. Barplots

#               3.1 for alkanes

#               3.2 for alkenes

#               3.3 for monomethylbranched CHC

#               3.4 for dimethylbranched CHC

#               3.5 for trimethylbranched CHC

#               3.6 for tetramethylbranched CHC

#               3.7 for dienes

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


# standardization with internal C12-standard to compensate deviations caused by injection volume

for (i in 1:nrow(dataset)) {
  
  dataset[i,]	<- dataset[i,]*22.5 / dataset[i,1]    # internal standard is here the first column 
  
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
                      levels = c("alkanes","alkenes","dienes","monomethyl-alkanes","dimethyl-alkanes","trimethyl-alkanes","tetramethyl-alkanes")) # sort group


p <- ggplot(dfall, aes(fill=group, y=val, x=caste)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 22) +  # white background 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(color = "black", size=14, face="bold")) +
  theme(legend.title =element_blank(), legend.text = element_text(color = "black", size = 16)) +
  scale_x_discrete(guide= guide_axis(n.dodge = 2)) +
  NULL +
  scale_fill_manual(values = c("goldenrod1", "darkorange2", "chartreuse2", "slategray1", "skyblue2", "steelblue3", "steelblue4")) +
  facet_grid(.~lineage, scales = "free")

p

# export graph
graph2ppt(p, file="Barplots-abs3_ST_20200826",width=15,height=11)
graph2png(p, file="Barplots-abs3_ST_20200826",width=15,height=11, dpi=600)
