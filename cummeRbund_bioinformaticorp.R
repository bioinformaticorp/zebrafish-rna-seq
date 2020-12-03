source('http://bioconductor.org/biocLite.R')
biocLite('cummeRbund')
library('cummeRbund')


setwd("~/curso/zebrafish/cuffdiff_out")
cuff<-readCufflinks()
cuff
genes<-genes(cuff)
disp<-dispersionPlot(genes)
disp # no funciona


#To assess the distributions of FPKM scores across samples, you can use the csDensity plot

#sin replicas
dens<-csDensity(genes(cuff))
dens

#con replicas
densRep<-csDensity(genes(cuff),replicates=T)
densRep

#Boxplots can be visualized using the csBoxplot method

#sin replicas
b<-csBoxplot(genes(cuff))
b

#con replicas
brep<-csBoxplot(genes(cuff),replicates=T)
brep

#A matrix of pairwise scatterplots can be drawn using the csScatterMatrix() method
s<-csScatterMatrix(genes(cuff))
s
#Individual Pairwise comparisons can be made by using csScatter
#You must specify the sample names to use for the x and y axes:

s<-csScatter(genes(cuff),"ZV9_2cells","ZV9_6h",smooth=T)
s

#Dendogram
#Sin repeticiones
dend<-csDendro(genes(cuff))

#Con repeticiones
dend.rep<-csDendro(genes(cuff),replicates=T)

#MvsA plots can be useful to determine any systematic bias that may be present between 
#conditions. The Cu Data method MAplot() can be used to examine these intensity vs 
#fold-change plots.  You must specify the sample names to use for the pairwise comparison 
#with x and y:
m<-MAplot(genes(cuff),"ZV9_2cells","ZV9_6h")
m
mCount<-MAplot(genes(cuff),"ZV9_2cells","ZV9_6h",useCount=T)
mCount


#Volcano plots are also available for the CuffData objects.
v<-csVolcanoMatrix(genes(cuff))
v

#For individual pairwise comparisons, you must specify the comparisons by sample name.
v<-csVolcano(genes(cuff),"ZV9_2cells","ZV9_6h")
v

##############################
#Cuffdiff  run information
##############################

#Run-level information such as run parameters,  and sample information can be accessed from a
#CuffSet object by using the runInfo and replicates methods:
runInfo(cuff)

#Cuando tienes replicas
replicates(cuff)



#######################################
#Creating Gene Sets
#######################################

geneIDs<-c("ENSDARG00000018303","ENSDARG00000040110","ENSDARG00000045067","ENSDARG00000053474","ENSDARG00000055868")
data(sampleData)
myGeneIds<-geneIDs
myGeneIds

myGenes<-getGenes(cuff,myGeneIds)
myGenes

#FPKM values for genes in gene set
head(fpkm(myGenes))

#Isoform-level FPKMs for gene set
head(fpkm(genes(myGenes)))

#Replicate FPKMs for TSS groups within gene set
head(repFpkm(TSS(myGenes)))


############################################
#Crear HEATMAPS
############################################

#Sin replicas
h<-csHeatmap(myGenes,cluster='both')
h

#Con replicas
h.rep<-csHeatmap(myGenes,cluster='both',replicates=T)
h.rep


#If you prefer barplots over heatmaps for genesets (although this is not necessarily 
#recommended for large gene sets).  You can use the expressionBarplot() method on a
#CuffFeatureSet or a CuffGeneSet object.
b<-expressionBarplot(myGenes)
b

#The csScatter() method can be used to produce scatter plot comparisons between any two 
#conditions.
s<-csScatter(myGenes,"ZV9_2cells","ZV9_6h",smooth=T)
s

#The volcano plot is a useful visualization to compare fold change between any two 
#conditions and signi cance (-log P-values).
v<-csVolcano(myGenes,"ZV9_2cells","ZV9_6h")
v

#Similar plots can be made for all sub-level features of a CuffGeneSet class by specifying 
#which slot you would like to plot (eg. isoforms(myGenes),TSS(myGenes),CDS(myGenes)).
ih<-csHeatmap(isoforms(myGenes),cluster='both',labRow=F)
ih
th<-csHeatmap(TSS(myGenes),cluster='both',labRow=F)
th

#Dendrograms can provide insight into the relationships between conditions for various 
#genesets (e.g.  significant genes used to draw relationships between conditions).  
#The method csDendro() can be used to plot a dendrogram based on Jensen-Shannon 
#distances between conditions for a given CuffFeatureSet or CuffGeneSet.
den<-csDendro(myGenes)

################################
#Significant genes alpha
###############################

#Graficar e identificar los genes con alpha 0.05
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)

mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')
mySigGeneIds
length(mySigGeneIds)

#Identificar los genes con alpha 0.05 por pares de condiciones
C1vsC2.sigIsoformIds<-getSig(cuff,x='ZV9_2cells',y='ZV9_6h',alpha=0.05,level='genes')
head(C1vsC2.sigIsoformIds)

#Heatmap entre condiciones
#Las similitudes entre las condiciones y / o las repeticiones pueden proporcionar 
#información útil sobre la relación entre varias agrupaciones de condiciones y pueden 
#ayudar a identificar las repeticiones atípicas que no se comportan como se esperaba.

#sin replica
myDistHeat<-csDistHeat(genes(cuff))

#con replica
myRepDistHeat<-csDistHeat(genes(cuff),replicates=T)



#####################################
#Demensionality reduction
####################################

#sin replica
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.MDS<-MDSplot(genes(cuff))

#Con replica
genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)
