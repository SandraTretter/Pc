# Sandra Tretter
# 19.08.2020


#### Pogonomyrmex californicus, CHC analysis, PCA

#---------------------------------------------------------------------

# content:  0. required packages/libraries

#           1. read in and transform the data

#           2. CHC classes

#           3. PCA

#               3.1 for all classes

#               3.2 for alkanes

#               3.3 for alkenes

#               3.4 for monomethylbranched CHCs

#               3.5 for dimethylbranched CHCs

#               3.6 for trimethylbranched CHCs

#               3.7 for tetramethylbranched CHCs

#               3.8 for dienes

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
#library(randomForest)
#library(reshape2)
#library(bipartite)
#library(dynRB)
#library(adegenet)
library(ggpubr)
library(export)


#---------------------------------------------------------------------

#### 1. read in and transform the data ####

setwd("C:/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten")

dataset <- read.csv("/Users/trett/Google Drive/AG Gadau/GC-MS/Pogonomyrmex/CHC analysis/Daten/BatchTable_CSV_ST_20200814.csv" ,header=T,dec=",",sep=";",check.names=FALSE,row.names = 1)


# Get rid of NAs

dataset[is.na(dataset)] = 0


# Create vector for "grouping" all samples 

vec <- vector(length=length(row.names(dataset))) 
for (i in 1:length(vec)) vec[i] = substr(row.names(dataset)[i],4,10) # Depends on sample names -> e.g. Pc_haplo_Q_1_01 has 15 characters -> 4,10 = Q_haplo


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

#### 2. CHC classes ####

nm <- colnames(dataset)

alkanes <- nm %in% grep("n-", nm, value = TRUE)
alkanes <- subset(dataset, select = alkanes)

alkenes <- nm %in% grep("alkene", nm, value = TRUE)
alkenes <- subset(dataset, select = alkenes)

me <- nm %in% grep("-Me", nm, value = TRUE)
me <- subset(dataset, select = me)

di <- nm %in% grep("-Di", nm, value = TRUE)
di <- subset(dataset, select = di)

tri <- nm %in% grep("-Tri", nm, value = TRUE)
tri <- subset(dataset, select = tri)

tetra <- nm %in% grep("-Tetra", nm, value = TRUE)
tetra <- subset(dataset, select = tetra)

diene <- nm %in% grep("diene", nm, value = TRUE)
diene <- subset(dataset, select = diene)

#-----------------------------------------------------------------------

#### 3. PCA ####

#-----------------------------------------------------------------------

#### 3.1 PCA for all classes ####

## with transformation

#clr-transformation
clr.output <- clr(dataset)
clr.output

dataset.clr.pca <- prcomp(clr.output, scale.=T, center=T)
summary(dataset.clr.pca)

all.clr <- ggbiplot(dataset.clr.pca, obs.scale=1, var.scale=1, choices=c(1,2), # choices -> which PCAs to plot
                var.axes=F, # remove arrows
                groups = vec, ellipse=TRUE) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  geom_point(aes(colour=vec, shape = vec, size = 3)) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) +
  ylim(-10,10) #+ 
  #geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2) # for sample names

all.clr 

graph2ppt(file="PCA2_ST_20200827.ppt",width = 13, height = 10)

## without transformation

dataset2.pca <- prcomp(dataset, scale.=T, center=T)
summary(dataset2.pca)

all <- ggbiplot(dataset2.pca, obs.scale=1, var.scale=1, choices=c(1,2), # choices -> which PCAs to plot
                 var.axes=F, # remove arrows
                 groups = vec, ellipse=TRUE) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  geom_point(aes(colour=vec, shape = vec, size = 3)) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) +
  ylim(-10,10) #+ 
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2) # for sample names

all

ggarrange(all.clr, all, labels = c("with clr", "without clr"), legend = "bottom", common.legend = TRUE)

#-----------------------------------------------------------------------

#### 3.2 PCA for alkanes #### 
clr.alkanes <- clr(alkanes)

alkanes.clr.pca <- prcomp(clr.alkanes, scale.=T, center=T)
summary(alkanes.clr.pca)

a.clr <- ggbiplot(alkanes.clr.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE, title("alkanes (decostand")) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

a.clr
graph2ppt(ae.clr, file="PCA-alkenes_ST_20200827.png",width=15,height=11)


alkanes.pca <- prcomp(alkanes, scale.=T, center=T)
summary(alkanes.pca)
a <- ggbiplot(alkanes.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE, title("alkanes (decostand")) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

a

ggarrange(a.clr, a, labels = c("alkanes with clr", "alkanes without clr"), legend = "bottom", common.legend = TRUE)

#-----------------------------------------------------------------------

#### 3.3 PCA for alkenes ####

clr.alkenes <- clr(alkenes)

alkenes.clr.pca <- prcomp(clr.alkenes, scale.=T, center=T)
summary(alkenes.clr.pca)

ae.clr <- ggbiplot(alkenes.clr.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) +
  ylim(c(-3,4)) + xlim(c(-7,4))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

ae.clr

alkenes.pca <- prcomp(alkenes, scale.=T, center=T)
summary(alkenes.pca)
ae <- ggbiplot(alkenes.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2))  +
  ylim(c(-3,4)) + xlim(c(-7,4))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

ae

ggarrange(ae.clr, ae, labels = c("alkenes with clr", "alkenes without clr"), legend = "bottom", common.legend = TRUE)
#-----------------------------------------------------------------------

#### 3.4 PCA for monomehtylbranched CHCs ####

clr.me <- clr(me)

me.clr.pca <- prcomp(clr.me, scale.=T, center=T)
summary(me.clr.pca)

m.clr <- ggbiplot(me.clr.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) +
  ylim(c(-7,5)) + xlim(c(-8,8))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

m.clr
graph2ppt(m.clr, file="PCA-mono_ST_20200827.png",width=15,height=11)

me.pca <- prcomp(me, scale.=T, center=T)
summary(me.pca)

m <- ggbiplot(me.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE, title("monomethylbranched (decostand")) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) +
  ylim(c(-7,5)) + xlim(c(-8,8))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

m

ggarrange(m.clr, m, labels = c("mono with clr", "mono without clr"), legend = "bottom", common.legend = TRUE)
#-----------------------------------------------------------------------

#### 3.5 PCA for dimethylbranched CHCs ####

clr.di <- clr(di)

di.clr.pca <- prcomp(clr.di, scale.=T, center=T)
summary(di.clr.pca)

d.clr <- ggbiplot(di.clr.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) +
  ylim(c(-7.5,6)) + xlim(c(-8,8.5))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

d.clr
graph2ppt(d.clr, file="PCA-di_ST_20200827.png",width=15,height=11)

di.pca <- prcomp(di, scale.=T, center=T)
summary(di.pca)

d <- ggbiplot(di.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE, title("monomethylbranched (decostand")) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) +
  ylim(c(-7.5,6)) + xlim(c(-8,8.5))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

d
ggarrange(d.clr, d, labels = c("di with clr", "di without clr"), legend = "bottom", common.legend = TRUE)
#-----------------------------------------------------------------------

#### 3.6 PCA for trimethylbranched CHCs ####

clr.tri <- clr(tri)

tri.clr.pca <- prcomp(clr.tri, scale.=T, center=T)
summary(tri.clr.pca)

t.clr <- ggbiplot(tri.clr.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) +
  ylim(c(-4,4)) + xlim(c(-6,6))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

t.clr
graph2ppt(t.clr, file="PCA-tri_ST_20200827.png",width=15,height=11)


tri.pca <- prcomp(tri, scale.=T, center=T)
summary(tri.pca)

t <- ggbiplot(tri.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE, title("monomethylbranched (decostand")) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) +
  ylim(c(-4,4)) + xlim(c(-6,6))
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

t
ggarrange(t.clr, t, labels = c("tri with clr", "tri without clr"), legend = "bottom", common.legend = TRUE)

#-----------------------------------------------------------------------

#### 3.7 PCA for tetramethylbranched CHCs ####

tetra.pca <- prcomp(tetra, scale.=T, center=T)
summary(tetra.pca)

tetra <- ggbiplot(tetra.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE, title("tetramethylbranched (decostand")) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) #+
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

tetra

#-----------------------------------------------------------------------

#### 3.8 PCA for dienes ####

diene.pca <- prcomp(diene, scale.=T, center=T)
summary(diene.pca)

diene <- ggbiplot(diene.pca, obs.scale=1, var.scale=1, choices=c(1,2), var.axes=F, groups = vec, ellipse=TRUE, title("dienes (decostand")) + 
  theme(legend.direction = 'vertical', legend.position = 'right', axis.title.x = element_text(color="black",size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme_classic() + 
  scale_color_manual(name=vec, values = c("skyblue1", "skyblue4", "violet", "violetred4", "orange", "orange4")) + 
  scale_shape_manual(name=vec, values=c(15, 0, 17, 2, 19, 1)) + 
  geom_point(aes(colour=vec, shape = vec, size = 2)) #+
#geom_text(aes(label=row.names(alkanes)), hjust=1, vjust=0.5, size=2)

diene

dev.off()
#-----------------------------------------------------------------------
