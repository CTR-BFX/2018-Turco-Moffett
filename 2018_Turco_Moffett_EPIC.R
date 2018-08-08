#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# EPIC Analysis        
# Data from GenomeScan EPIC Array (organoids)
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2018-Turco-Moffett
#
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics Facility
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

# Load in the required R packages
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("minfi","ChAMPdata","Illumina450ProbeVariants.db","sva","IlluminaHumanMethylation450kmanifest","limma","RPMM","DNAcopy","preprocessCore","impute","marray","wateRmelon","goseq","plyr","GenomicRanges","RefFreeEWAS","qvalue","isva","doParallel","bumphunter","quadprog","shiny","shinythemes","plotly","RColorBrewer","DMRcate","dendextend","IlluminaHumanMethylationEPICmanifest","FEM","matrixStats","missMethyl","combinat"))
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library("ggplot2")
library('minfi')
library('lumi')
#library('wateRmelon')
library('methylumi')
library('limma')
library('VennDiagram');
#library('methyAnalysis')
require('FDb.InfiniumMethylation.hg19')
require('IlluminaHumanMethylationEPICmanifest')
#require('IlluminaHumanMethylation450k.db')
library('ChAMP')
library("data.table")
library("cowplot")
library("reshape")
library("ggrepel")
library("ggdendro")
library("useful")
library("MASS")
library("viridis")
library("scales")
library("biomaRt")
require("plyr")
library("ggalt")


#ensembl    <-  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), mart = ensembl)

#grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
#listAttributes(grch37)
# http://zwdzwd.github.io/InfiniumAnnotation
#test <- read.table(gzfile("EPIC.hg38.manifest.gencode.v22.tsv.gz"), sep="\t", header=TRUE)
#head(test, 100)

#library("IlluminaHumanMethylationEPICmanifest")
#library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

#
# Load in the EPIC probe manifest
#
data(probe.features.epic)
str(probe.features)
probe.features$probe_id <- rownames(probe.features)
head(probe.features)


probe.locations.epic <- probe.features[, c("CHR","MAPINFO", "Strand")]
probe.locations.epic$probe_id <- rownames(probe.locations.epic)
head(probe.locations.epic)


#------------------------------------------------------------------------------
# Set the directories and baseDir for the 450K experiment
setwd("/Users/rhamilto/Documents/CTR-Data/CTR_EPIC/")
baseDir <- getwd()
baseDirIDATs <- sprintf("%s/IDATs",baseDir)

Project <- "CTR_EPIC.all"
#------------------------------------------------------------------------------

IDATs.MB             <- "/Users/rhamilto/Documents/CTR-Data/CTR_EPIC/E-GEOD-66210"
myLoad.MB            <- champ.load(directory = IDATs.MB, methValue='B', filterBeads=TRUE, filterXY=TRUE, detPcut = 0.01, arraytype="450K")
head(myLoad.MB$pd)

sampleTable.MB.champ <- as.data.frame(myLoad.MB$pd)
myQC.MB              <- champ.QC(beta=myLoad.MB$beta, pheno=myLoad.MB$pd$Sample_Group, resultsDir="./CHAMP_QCimages_MB/")

myNorm.MB            <- champ.norm(beta=myLoad.MB$beta, rgSet=myLoad.MB$rgSet, mset=myLoad.MB$mset, resultsDir="./CHAMP_Normalization_MB/")


#------------------------------------------------------------------------------

myLoad <- champ.load(directory = baseDirIDATs, methValue='B', filterBeads=TRUE, filterXY=TRUE, detPcut = 0.01, arraytype="EPIC")
head(myLoad$pd)

sampleTable.champ <- as.data.frame(myLoad$pd)

myQC   <- champ.QC()

myNorm <- champ.norm()


myLoad <- ""



customPCADendro <- function(ProjectTitle, myNorm, TOPNUM, sampleTable.champ) {
  
  myNorm.df <- as.data.frame(myNorm)
  rv        <- rowVars(as.matrix(myNorm.df))
  select    <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca       <- prcomp(t(myNorm.df[select, ]))
  pc1var    <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var    <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab    <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab    <- paste0("PC2 (",as.character(pc2var),"%)")
  scores    <- data.frame(sampleName=sampleTable.champ$Sample_Name, pca$x, Sample_Group=sampleTable.champ$Sample_Group, Pool_ID=sampleTable.champ$Pool_ID)

  scores$sampleName <- gsub("_BS", "", scores$sampleName)
  scores$sampleName <- gsub("_oxBS", "", scores$sampleName)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, col = factor(Sample_Group), shape=factor(Pool_ID)) ) +
             geom_point(size = 3, alpha=0.75 ) + 
             geom_text_repel(aes(label=sampleName), show.legend = FALSE) +
             xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
             ggtitle(paste0(ProjectTitle, " PCA Top ", TOPNUM, " MV")) +
             scale_colour_manual(name="Sample Group", values = c("control"="red", "organoid"="green", "preflow"="purple")) +
             scale_shape_manual(name="Treatment", values = c(19,17)) +
             theme(text = element_text(size=12)) 
  
  png(paste0(Project, ".QC.PCA.", TOPNUM, '.png'), units="cm", width=15, height=15, res=250)
  par(bg=NA)
  print(plt.pca)
  dev.off()
  
  dd.row     <- as.dendrogram(hclust(dist(t( myNorm.df[select, ]))))
  ddata      <- dendro_data(dd.row)
  labs       <- label(ddata)
  labs$group <- labs$label
  labs$group <- gsub("pre.*",  "preflow", labs$group)
  labs$group <- gsub("org.*",  "organoid", labs$group)
  labs$group <- gsub("cont.*", "control", labs$group)
  

  plt.dendro <- ggplot(segment(ddata)) +
                geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                geom_text(data=label(ddata),aes(label=label, x=x, y=-0.25, angle=-90, colour=labs$group), hjust=0) +
                scale_colour_manual(name="Sample Group", values = c("control"="red", "organoid"="green", "preflow"="purple"), guide=FALSE) +
                ylim(-3, max(ddata$segments$y)) + ylab("distance") + xlab("") +
                ggtitle(paste0(ProjectTitle, " Dendrogram Top ", TOPNUM, " MV")) +
                theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank() )
    
  png(paste0(Project, ".QC.dendrogram.", TOPNUM, '.png'), units="cm", width=15, height=15, res=250)
  par(bg=NA)
  print(plt.dendro)
  dev.off()
  
  myNorm.df <- ""
  
  return(list(plt.pca, plt.dendro))
}
plt.QC <- customPCADendro(Project, myNorm, 500, sampleTable.champ)
plt.QC[[1]]
plt.QC[[2]]


functionPlotCorrelation_BSoxBS <- function(myNorm, idx1, idx2, threshold, sampleTable.champ) {
  
  myNorm.df           <- as.data.frame(myNorm)
  colnames(myNorm.df) <- sampleTable.champ$Sample_Name
  x                   <- myNorm.df[,colnames(myNorm.df)[idx2]]
  y                   <- myNorm.df[,colnames(myNorm.df)[idx1]]
  xy                  <- data.frame(x,y)
  x <- ""
  y <- ""
  xy$diff             <- (xy$y - xy$x)
  xt                  <- colnames(myNorm.df)[idx2]
  yt                  <- colnames(myNorm.df)[idx1]
  r2                  <- paste("R^2 == ", as.character(format(summary(lm(y ~ x, xy))$r.squared, digits=5)))

  myNorm.df <- ""
  
  # http://slowkow.com/notes/ggplot2-color-by-density/
  get_density <- function(x, y, n = 100) {
                          dens <- MASS::kde2d(x = x, y = y, n = n)
                          ix   <- findInterval(x, dens$x)
                          iy   <- findInterval(y, dens$y)
                          ii   <- cbind(ix, iy)
                          return(dens$z[ii]) }
  
  xy$density <- get_density(xy$x, xy$y)
  
  plt <- ggplot(xy, aes(x=x, y=y)) + 
         geom_point(alpha=0.25, size=0.1, aes(x=x, y=y, colour=density)) +
         scale_color_viridis() +
         geom_abline(intercept = threshold, slope = 1, linetype="dashed") +
         geom_abline(intercept = 0, slope = 1, linetype="dashed") +
         geom_abline(intercept = -(threshold), slope = 1, linetype="dashed") +
         annotate("text", label = r2, x = 0.2, y = 1, size = 2, colour = "blue", parse=TRUE) +
         geom_smooth(method=lm, colour = "blue", linetype = "solid", size=0.5, alpha=0.95) +
         coord_fixed() +
         ylab("BS (mC+hmC)") + xlab("oxBS (mC)") +
         ggtitle( bquote(.(yt)~vs~.(xt)~"[threshold"~Delta*beta ==.(threshold)*"]") ) +
         theme(text=element_text(size=6, family="sans"), plot.title=element_text(size=6, family="sans"),
               axis.text.x=element_text(size=6, family="sans"), axis.text.y=element_text(size=6, family="sans"),
               legend.position="none", 
               legend.text=element_text(size=4), legend.title=element_text(size=6), legend.key.height = unit(0.2, "cm"))
  
  xy        <- ""
  
  png(paste0("SamplePlots/", Project, ".BSoxBS.Correlation.", yt, ".", xt, '.png'), units="cm", width=7.5, height=7.5, res=250)
  par(bg=NA)
  print(plt)
  dev.off()
  
  return(plt)
}
functionPlotCorrelation_BSoxBS(myNorm, 3, 4, 0.2, sampleTable.champ)

 

functionPlotBeta_BSoxBS <- function(myNorm, idx1, idx2, sampleTable.champ) {
  
  myNorm.df           <- as.data.frame(myNorm)
  colnames(myNorm.df) <- sampleTable.champ$Sample_Name
  x                   <- myNorm.df[,colnames(myNorm.df)[idx1]]
  y                   <- myNorm.df[,colnames(myNorm.df)[idx2]]
  xt                  <- colnames(myNorm.df)[idx1]
  yt                  <- colnames(myNorm.df)[idx2]
  xy                  <- data.frame(x,y)
  colnames(xy)        <- c("BS", "oxBS")
  
  suppressMessages({
  xy.mlt              <- melt(xy) #, c("BS", "oxBS") )
  })
  colnames(xy.mlt)    <- c("Treatment", "Beta")

  delta.lab <- as.character(bquote(Delta))
  beta.lab  <- as.character(bquote(beta))
  
  plt <- ggplot(data=xy.mlt) + 
         geom_line( aes(Beta, colour=Treatment, fill=Treatment), stat="density", alpha=0.9 ) +
         scale_colour_manual(name="", values=c("BS"="black", "oxBS"="purple")) +
         ggtitle(paste0(yt, " vs ", xt)) +
         theme(text=element_text(size=8, family="sans"), plot.title=element_text(size=6, family="sans"),
               legend.position=c(0.5, 0.85),
               axis.text.x=element_text(size=8, family="sans"), axis.text.y=element_text(size=8, family="sans"))
         guides(colour = guide_legend(override.aes = list(size=5))) 
  
  png(paste0("SamplePlots/", Project, ".BSoxBS.Beta.", yt, ".", xt, '.png'), units="cm", width=7.5, height=7.5, res=250)
  par(bg=NA)
  print(plt)
  dev.off()
  
  myNorm.df <- ""
  
  return(plt)
}
functionPlotBeta_BSoxBS(myNorm, 3, 4, sampleTable.champ)


functionPlotBetaDistribution <- function(myNorm, idx1, idx2, threshold, limit, sampleTable.champ) {
  
  myNorm.df           <- as.data.frame(myNorm)
  colnames(myNorm.df) <- sampleTable.champ$Sample_Name
  x                   <- myNorm.df[,colnames(myNorm.df)[idx2]]
  y                   <- myNorm.df[,colnames(myNorm.df)[idx1]]
  xy                  <- data.frame(x,y)
  xy$diff             <- (xy$y - xy$x)
  xt                  <- colnames(myNorm.df)[idx1]
  yt                  <- colnames(myNorm.df)[idx2]

  xy.inset            <- subset(xy, abs(xy$diff) > limit)
  
  plt <- ggplot(xy, aes(x=diff)) + 
         geom_histogram(aes(y=..density..), binwidth=.01, colour="black", fill="white", size=0.25) +
         geom_density(alpha=.2, fill="#FF6666") +
         geom_vline(xintercept = -(threshold), linetype="dashed") +
         geom_vline(xintercept = 0) +
         geom_vline(xintercept = threshold, linetype="dashed") +
         xlab(bquote(Delta*beta)) +
         ggtitle( bquote(atop(.(yt)~vs~.(xt),"[threshold"~Delta*beta ==.(threshold)*"]")) ) +
         theme(text=element_text(size=6, family="sans"), plot.title=element_text(size=6, family="sans"),
               axis.text.x=element_text(size=6, family="sans"), axis.text.y=element_text(size=6, family="sans"),
               legend.position="none")

  plt.inset <- ggplot(xy.inset, aes(x=diff)) + 
               geom_histogram(aes(y=..density..), binwidth=.01, colour="black", fill="white", size=0.25) +
               geom_density(alpha=.2, fill="#FF6666") +
               geom_vline(xintercept = -(threshold), linetype="dashed") +
               geom_vline(xintercept = 0) +
               geom_vline(xintercept = threshold, linetype="dashed") +
               xlab(bquote(Delta*beta)) + 
               #ylim(0,1.0) + 
               ggtitle( bquote(atop("Inset"~.(yt)~vs~.(xt),"[threshold"~Delta*beta ==.(threshold)*"]") )) +
               theme(text=element_text(size=6, family="sans"), plot.title=element_text(size=6, family="sans"),
               axis.text.x=element_text(size=6, family="sans"), axis.text.y=element_text(size=6, family="sans"),
               legend.position="none")

  cow.plt <- plot_grid(plt, plt.inset, labels=c("A", "B"), rel_widths = c(1,1))
  
  myNorm.df <- ""
  xy        <- ""
  xy.inset  <- ""
  
  png(paste0("SamplePlots/", Project, ".BSoxBS.DeltaBeta.", yt, ".", xt, '.png'), units="cm", width=10, height=7.5, res=250)
  par(bg=NA)
  print(cow.plt)
  dev.off()
  
  return(cow.plt)
  }
functionPlotBetaDistribution(myNorm, 3, 4, 0.2, 0.2, sampleTable.champ)


functionPlotBetaPercentage <- function(myNorm, idx1, sampleTable.champ) {
  myNorm.df           <- as.data.frame(myNorm)
  colnames(myNorm.df) <- sampleTable.champ$Sample_Name
  x                   <-  myNorm.df[,colnames(myNorm.df)[idx1]]
  x                   <- as.data.frame(x)
  xt                  <- colnames(myNorm.df)[idx1]

  plt <- ggplot(x,aes(x=x)) +
         geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=.2, boundary=0, 
                        colour="black", fill="grey", alpha=0.75, size=0.1) +
         scale_y_continuous(labels = percent, breaks=seq(0,1,0.05)) +
         scale_x_continuous(breaks=seq(0,1,0.2)) +
         xlab(bquote(beta~"(methylation)")) + 
         ylab("Percentage") +
         ggtitle(paste0(xt)) +
         theme(text=element_text(size=6, family="sans"), plot.title=element_text(size=6, family="sans"),
               axis.text.x=element_text(size=6, family="sans"), axis.text.y=element_text(size=6, family="sans"),
               legend.position="none")
 
  png(paste0("SamplePlots/", Project, ".BSoxBS.PercentBeta.", xt, '.png'), 
      units="cm", width=7.5, height=7.5, res=250)
  par(bg=NA)
  print(plt)
  dev.off()
return(plt)
}
functionPlotBetaPercentage(myNorm, 3, sampleTable.champ)




counter = 1
corr.plot.list <- list() 
beta.plot.list <- list() 
for (i in seq(1, (nrow(sampleTable.champ)), 2))
{
  message(paste0(counter, ": Making BS/oxBS correlation plot for samples ", i, " and ", i+1))
  corr.plt <- functionPlotCorrelation_BSoxBS(myNorm, i, i+1, 0.2, sampleTable.champ)

  message(paste0(counter, ": Making BS/oxBS beta plot for samples ", i, " and ", i+1))
  beta.plt <- functionPlotBeta_BSoxBS(myNorm, i, i+1, sampleTable.champ)

  message(paste0(counter, ": Making BS/oxBS Delta beta plot for samples ", i, " and ", i+1))
  delta.beta.plt <- functionPlotBetaDistribution(myNorm, i, i+1, 0.2, 0.15, sampleTable.champ)
  
  message(paste0(counter, ": Making BS/oxBS Percent beta plot for samples ", i, " and ", i+1))
  perc.delta.plt <- functionPlotBetaPercentage(myNorm, i,   sampleTable.champ)
  perc.delta.plt <- functionPlotBetaPercentage(myNorm, i+1, sampleTable.champ)
  
  counter <- counter + 1
}


#
# Look at the genomic positions
#
threshold <- 0.8


data(probe.features.epic)
str(probe.features)
probe.features$probe_id <- rownames(probe.features)
head(probe.features)
nrow(probe.features)


#
# Annotate the EPIC organoid samples
#
nrow(myNorm)
myNorm.probes.df                <- as.data.frame(myNorm)
myNorm.probes.df$probe_id       <- rownames(myNorm.probes.df)
myNorm.probes.ann               <- merge(myNorm.probes.df, probe.features, by="probe_id")
myNorm.probes.ann$Strand        <- gsub("F", "1", myNorm.probes.ann$Strand)
myNorm.probes.ann$Strand        <- gsub("R", "-1", myNorm.probes.ann$Strand)

myNorm.probes.ann$promoter      <- myNorm.probes.ann$feature
myNorm.probes.ann$promoter      <- gsub("TSS1500", "promoter", myNorm.probes.ann$promoter)
myNorm.probes.ann$promoter      <- gsub("TSS200", "promoter", myNorm.probes.ann$promoter)

myNorm.probes.ann$org_meanBeta  <- rowMeans( myNorm.probes.ann[,grep("org",colnames(myNorm.probes.ann))] )
myNorm.probes.ann$pre_meanBeta  <- rowMeans( myNorm.probes.ann[,grep("pre",colnames(myNorm.probes.ann))] )
myNorm.probes.ann$cntl_meanBeta <- rowMeans( myNorm.probes.ann[,grep("control",colnames(myNorm.probes.ann))] )
head(myNorm.probes.ann)
nrow(myNorm.probes.ann)


#
# Annotate the 450K maternal Blood samples
#
nrow(myNorm.MB)
myNorm.MB.probes.df                <- as.data.frame(myNorm.MB)
myNorm.MB.probes.df$probe_id       <- rownames(myNorm.MB.probes.df)
myNorm.MB.probes.ann               <- merge(myNorm.MB.probes.df, probe.features, by="probe_id")
myNorm.MB.probes.ann$Strand        <- gsub("F", "1", myNorm.MB.probes.ann$Strand)
myNorm.MB.probes.ann$Strand        <- gsub("R", "-1", myNorm.MB.probes.ann$Strand)

myNorm.MB.probes.ann$promoter      <- myNorm.MB.probes.ann$feature
myNorm.MB.probes.ann$promoter      <- gsub("TSS1500", "promoter", myNorm.MB.probes.ann$promoter)
myNorm.MB.probes.ann$promoter      <- gsub("TSS200", "promoter", myNorm.MB.probes.ann$promoter)

myNorm.MB.probes.ann$mb_meanBeta  <- rowMeans( myNorm.MB.probes.ann[,grep("mb_",colnames(myNorm.MB.probes.ann))] )

head(myNorm.MB.probes.ann)
nrow(myNorm.MB.probes.ann)



functionGetPerSampleFeature <- function(myNorm, Feature, SubFeature, sampleGroup) {
  
  myNorm           <- melt(myNorm[, c(Feature, colnames(myNorm)[grep("meanBeta",colnames(myNorm))] )])
  colnames(myNorm) <- c("Feature","variable","value")
  myNorm           <- subset(myNorm, Feature==SubFeature)

  print( head(myNorm))
  
  return(myNorm) 
}


cgi.island.org  <- functionGetPerSampleFeature(myNorm.probes.ann, "cgi", "island", "org_meanBeta")  
cgi.island.mb   <- functionGetPerSampleFeature(myNorm.MB.probes.ann, "cgi", "island", "mb_meanBeta")  


feature.body.org  <- functionGetPerSampleFeature(myNorm.probes.ann, "feature", "Body", "org_meanBeta")
feature.body.mb  <- functionGetPerSampleFeature(myNorm.MB.probes.ann, "feature", "Body", "mb_meanBeta")


promoter.org   <- functionGetPerSampleFeature(myNorm.probes.ann, "promoter", "promoter", "org_meanBeta")  
promoter.mb   <- functionGetPerSampleFeature(myNorm.MB.probes.ann, "promoter", "promoter", "mb_meanBeta")  



#
# L1 Overlays
#
mb_1_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_1_BS.intersection.L1.bed"))
mb_2_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_2_BS.intersection.L1.bed"))
mb_3_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_3_BS.intersection.L1.bed"))
mb_4_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_4_BS.intersection.L1.bed"))
mb_5_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_5_BS.intersection.L1.bed"))
L1.mb                     <- cbind(mb_1_BS.intersection.L1,mb_2_BS.intersection.L1$V4, mb_3_BS.intersection.L1$V4,mb_4_BS.intersection.L1$V4, mb_4_BS.intersection.L1$V4)
L1.mb$value               <- rowMeans(L1.mb[, c(4,5,6,7,8)])
L1.mb$Feature             <- "L1"
L1.mb$variable            <- "mb_meanBeta"
L1.mb                     <- L1.mb[, c("Feature", "variable", "value")]
head(L1.mb)

control_BS.intersection.L1 <- read.table(paste0(baseDir, "/FeatureBedFiles/", "control_BS.intersection.L1.bed"))
L1.control                 <- control_BS.intersection.L1
L1.control$value           <- L1.control[, c("V4")]
L1.control$Feature         <- "L1"
L1.control$variable        <- "cntl_meanBeta"
L1.control                 <- L1.control[, c("Feature", "variable", "value")]
head(L1.control)

org_1_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_1_BS.intersection.L1.bed"))
org_3_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_3_BS.intersection.L1.bed"))
org_4_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_4_BS.intersection.L1.bed"))
org_5_BS.intersection.L1   <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_5_BS.intersection.L1.bed"))
L1.org                     <- cbind(org_1_BS.intersection.L1,org_3_BS.intersection.L1$V4,org_4_BS.intersection.L1$V4,org_5_BS.intersection.L1$V4)
L1.org$value               <- rowMeans(L1.org[, c(4,5,6,7)])
L1.org$Feature             <- "L1"
L1.org$variable            <- "org_meanBeta"
L1.org                     <- L1.org[, c("Feature", "variable", "value")]
head(L1.org)

pre_63_BS.intersection.L1 <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_63_BS.intersection.L1.bed"))
pre_64_BS.intersection.L1 <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_64_BS.intersection.L1.bed"))
pre_65_BS.intersection.L1 <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_65_BS.intersection.L1.bed"))
pre_66_BS.intersection.L1 <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_66_BS.intersection.L1.bed"))
L1.pre                    <- cbind(pre_63_BS.intersection.L1,pre_64_BS.intersection.L1$V4,pre_65_BS.intersection.L1$V4,pre_66_BS.intersection.L1$V4)
L1.pre$value              <- rowMeans(L1.pre[, c(4,5,6,7)])
L1.pre$Feature            <- "L1"
L1.pre$variable           <- "pre_meanBeta"
L1.pre                    <- L1.pre[, c("Feature", "variable", "value")]
head(L1.pre)





feature.summary.tbl          <- rbind(cgi.island.org,   cgi.island.mb,   
                                      promoter.org,     promoter.mb,     
                                      feature.body.org, feature.body.mb, 
                                      L1.org,   L1.pre,   L1.control, L1.mb )
head(feature.summary.tbl)
unique(feature.summary.tbl$variable)
unique(feature.summary.tbl$Feature)


pdf(paste(Project, "-GenomicFeatureSummary_Org_vs_Pre_Version4.pdf", sep=""), width=6,height=5)
par(bg=NA)
ggplot(data=feature.summary.tbl, aes(x=Feature, y=value)) +
  geom_boxplot(aes(fill=variable), outlier.alpha = 0) +
  scale_fill_manual(name="", values = c("org_meanBeta"="lightpink1", "pre_meanBeta"="dodgerblue", "cntl_meanBeta"="grey", "mb_meanBeta"="red"), 
                    limits=c("org_meanBeta","pre_meanBeta","cntl_meanBeta","mb_meanBeta"), labels=c("organoid","placenta(pre-flow)","control/brain", "maternal blood")) +
  scale_x_discrete(name="", labels=c("CpG Islands", "Promoters", "Gene Bodies", "LINE1")) +
  ylab( bquote("Methylation ("*beta~"value)")) +
  ylim(0,1) +
  theme(text=element_text(size=12,  family="sans"), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="top")
dev.off()




#
# Specific Gene
#

functionGetGeneMethylationStats <- function(GENE) {
 # GENE <- "ELF5"  
  
  mb_1_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_1_BS.intersect.", GENE, ".bed"))
  mb_2_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_2_BS.intersect.", GENE, ".bed"))
  mb_3_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_3_BS.intersect.", GENE, ".bed"))
  mb_4_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_4_BS.intersect.", GENE, ".bed"))
  mb_5_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "mb_5_BS.intersect.", GENE, ".bed"))
  Gene.mb                     <- cbind(mb_1_BS.intersection,mb_2_BS.intersection$V4, mb_3_BS.intersection$V4,mb_4_BS.intersection$V4,mb_5_BS.intersection$V4)
  Gene.mb$mb_meanBeta         <- rowMeans(Gene.mb[, c(4,5,6,7,8)])
  Gene.mb.summary             <- summary(Gene.mb$mb_meanBeta)
  Gene.mb.summary["Std.Dev."] <- round(sd(Gene.mb$mb_meanBeta),7)
  Gene.mb.summary             <- data.frame(sample="mb_meanBeta", feature=GENE, 
                                             ymin=Gene.mb.summary[[1]], lower=Gene.mb.summary[[2]], 
                                             middle=Gene.mb.summary[[3]], 
                                             upper=Gene.mb.summary[[5]], ymax=Gene.mb.summary[[6]], 
                                             sd=Gene.mb.summary[[7]], mean=Gene.mb.summary[[4]])
  Gene.mb.summary$ymin[Gene.mb.summary$ymin < 0] <- 0
  Gene.mb.summary$ymax[Gene.mb.summary$ymax > 1] <- 1  
 
   
org_1_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_1_BS.intersect.", GENE, ".bed"))
org_3_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_3_BS.intersect.", GENE, ".bed"))
org_4_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_4_BS.intersect.", GENE, ".bed"))
org_5_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "org_5_BS.intersect.", GENE, ".bed"))
Gene.org                     <- cbind(org_1_BS.intersection,org_3_BS.intersection$V4,org_4_BS.intersection$V4,org_5_BS.intersection$V4)
Gene.org$org_meanBeta        <- rowMeans(Gene.org[, c(4,5,6,7)])

Gene.org.summary             <- summary(Gene.org$org_meanBeta)
Gene.org.summary["Std.Dev."] <- round(sd(Gene.org$org_meanBeta),7)
Gene.org.summary             <- data.frame(sample="org_meanBeta", feature=GENE, 
                                           ymin=Gene.org.summary[[1]], lower=Gene.org.summary[[2]], 
                                           middle=Gene.org.summary[[3]], 
                                           upper=Gene.org.summary[[5]], ymax=Gene.org.summary[[6]], 
                                           sd=Gene.org.summary[[7]], mean=Gene.org.summary[[4]])
Gene.org.summary$ymin[Gene.org.summary$ymin < 0] <- 0
Gene.org.summary$ymax[Gene.org.summary$ymax > 1] <- 1


pre_63_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_63_BS.intersect.", GENE, ".bed"))
pre_64_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_64_BS.intersect.", GENE, ".bed"))
pre_65_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_65_BS.intersect.", GENE, ".bed"))
pre_66_BS.intersection        <- read.table(paste0(baseDir, "/FeatureBedFiles/", "pre_66_BS.intersect.", GENE, ".bed"))
Gene.pre                     <- cbind(pre_63_BS.intersection,pre_64_BS.intersection$V4,pre_65_BS.intersection$V4,pre_66_BS.intersection$V4)
Gene.pre$pre_meanBeta        <- rowMeans(Gene.pre[, c(4,5,6,7)])
Gene.pre.summary             <- summary(Gene.pre$pre_meanBeta)
Gene.pre.summary["Std.Dev."] <- round(sd(Gene.pre$pre_meanBeta),7)
Gene.pre.summary             <- data.frame(sample="pre_meanBeta", feature=GENE, 
                                           ymin=Gene.pre.summary[[1]], lower=Gene.pre.summary[[2]], 
                                           middle=Gene.pre.summary[[3]], 
                                           upper=Gene.pre.summary[[5]], ymax=Gene.pre.summary[[6]], 
                                           sd=Gene.pre.summary[[7]], mean=Gene.pre.summary[[4]])
Gene.pre.summary$ymin[Gene.pre.summary$ymin < 0] <- 0
Gene.pre.summary$ymax[Gene.pre.summary$ymax > 1] <- 1


Gene.data           <- cbind(Gene.org, Gene.pre$pre_meanBeta) 
Gene.data$Gene      <- GENE
Gene.data           <- Gene.data[, c("Gene", "V2", "org_meanBeta", "Gene.pre$pre_meanBeta")]
colnames(Gene.data) <- c("Gene", "V2", "org_meanBeta", "pre_meanBeta")

Gene.mb$Gene        <- GENE
Gene.mb             <- Gene.mb[, c("Gene", "V2", "mb_meanBeta")]

Gene.data           <- merge(Gene.data,Gene.mb, all=TRUE) # For EPIC/450K probes

Gene.data           <- Gene.data[, c("Gene", "V2", "org_meanBeta", "pre_meanBeta", "mb_meanBeta")]
colnames(Gene.data) <- c("Gene", "CpG", "Org", "PL", "MB")

Gene.data           <- melt(Gene.data, id.vars=c("CpG", "Gene"))

return(list(Gene.mb,Gene.org, Gene.pre, Gene.data))
}


selected.gene.ELF5    <- functionGetGeneMethylationStats("ELF5")
selected.gene.EZR     <- functionGetGeneMethylationStats("EZR")
selected.gene.HAND1   <- functionGetGeneMethylationStats("HAND1")
selected.gene.LASP1   <- functionGetGeneMethylationStats("LASP1")
selected.gene.MAP3K8  <- functionGetGeneMethylationStats("MAP3K8")
selected.gene.PLET1   <- functionGetGeneMethylationStats("PLET1")
selected.gene.RIN3    <- functionGetGeneMethylationStats("RIN3")
selected.gene.SH2D3C  <- functionGetGeneMethylationStats("SH2D3C")
selected.gene.TEAD4   <- functionGetGeneMethylationStats("TEAD4")
selected.gene.TINAGL1 <- functionGetGeneMethylationStats("TINAGL1")

GeneData.tbl <- rbind(selected.gene.ELF5[[4]],  selected.gene.EZR[[4]],    selected.gene.HAND1[[4]], 
                      selected.gene.LASP1[[4]], selected.gene.MAP3K8[[4]], selected.gene.PLET1[[4]], 
                      selected.gene.RIN3[[4]],  selected.gene.SH2D3C[[4]], selected.gene.TEAD4[[4]],
                      selected.gene.TINAGL1[[4]])
head(GeneData.tbl)

pdf(paste(Project, "-SelectedGenes_Org_vs_Pre_Version2.pdf", sep=""), width=7,height=5)
par(bg=NA)
ggplot(data=GeneData.tbl, aes(x=variable, y=value, group=Gene, colour=variable, fill=variable)) +
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = .4, fill=NA) +
  geom_jitter(pch = 21, width=0.1, size=0.5) +
  scale_colour_manual(name="", values = c("Org"="lightpink1", "PL"="dodgerblue", "MB"="red")) +
#  ylim(0,1) + 
  xlab("") +
  ylab( bquote("Methylation ("*beta~"value)")) +
  facet_wrap(~ Gene, nrow = 2) +
  ggtitle("Gene Body Methylation (Individual CpGs)") +
  theme(text=element_text(size=12,  family="sans"), #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(color="red", size=0.25, fill=NA),
        legend.position="none")
dev.off()




MYLOAD  <- myLoad.MB$pd
MYNORM  <- myNorm.MB.probes.ann



for (i in seq(1, (nrow(MYLOAD)), 1))
{
  sampleName <- colnames(MYNORM)[i+1] 
  print( paste0("Generating BED/BedGraph for ", sampleName))
  
  BED <- MYNORM[, c(nrow(MYLOAD)+2, nrow(MYLOAD)+3, nrow(MYLOAD)+3, 1, i+1, nrow(MYLOAD)+4)]

  #chromA  chromStartA  chromEndA  dataValueA
  BG  <- MYNORM[, c(nrow(MYLOAD)+2, nrow(MYLOAD)+3, nrow(MYLOAD)+3, i+1) ]
  #BG  <- head(BG, 10)
  
  write.table(BED, file=paste0("Methylation_Files/", sampleName, ".bed"),      col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  write.table(BG,  file=paste0("Methyaltion_Files/", sampleName, ".bedgraph"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  
  fConn   <- file(paste0("Methyaltion_Files/", sampleName, ".bedgraph"), 'r+') 
  Lines   <- readLines(fConn) 
  topline <- paste0("track type=bedGraph name=", sampleName, " description=bedGraph visibility=full color=200,100,0 altColor=0,100,200 priority=20")
  writeLines(c(topline, Lines), con = fConn) 
  close(fConn) 
}
  
  
  
  
  
cg.dot.tbl                  <- read.table("ELF5_Promoter.dot", header = TRUE, sep=",")
cg.dot.tbl[cg.dot.tbl == 0] <- NA
cg.dot.tbl.mlt              <- melt(cg.dot.tbl, id=c("CG"))
cg.dot.tbl.mlt$variable     <- gsub(".srt.bedgraph", "", cg.dot.tbl.mlt$variable)
head(cg.dot.tbl.mlt, 10)

pdf(paste(Project, "-ELF5_Promoter_Meth.pdf", sep=""), width=7,height=3.0)
par(bg=NA)
ggplot(cg.dot.tbl.mlt, aes(y = variable, x = CG)) +
  geom_point(aes(fill=value), size=6, colour="black", shape=21) +
  scale_fill_gradient(name=bquote(atop("Methylation",beta~"value")), low = "white", high = "black", limits=c(0,1), space = "Lab", na.value = "pink", guide = "colourbar", aesthetics = "fill") +
  scale_x_reverse(name="CG Position", breaks=seq(1,11,1), labels=seq(11,1,-1)) +
  xlab("CG Position") + ylab("") + 
  ggtitle("ELF5 Promoter")
dev.off()




#------------------------------------------------------------------------------
# END OF SCRIPT
#------------------------------------------------------------------------------