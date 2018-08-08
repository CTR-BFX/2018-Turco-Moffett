#!/usr/local/bin/Rscript
rm(list=ls())  # Read more at: http://scq.io/AVkLvbhX#gs.WJrwyFg

#
# CTR_myt25_0002 ::: Microarray Analysis ::: Analysis on normalised data  
# Copyright Malwina Prater (mn367@cam.ac.uk)
#


# The aim of this microarray experiment is to:
# (1) find out how closely they represent the trophoblast cells, in vivo
# So far, I have done a simple comparison of Stroma vs TOrg and many interesting syncytial and trophoblast genes (as well as imprinted genes) that we know, came up
# (2) learn more about the genes expressed and how they describe the trophoblast cell identity
# Many interesting LNC RNAs and transcription factors came up as well. I was wondering whether we could use this information to learn more about the signalling pathways that are overrepresented in trophoblast organoids.
# It would be great if we could think and discuss about how to mine this data further and to find the best way to represent the analysis for the paper.

library("data.table")
library("pheatmap")
library("grid")
library("biomaRt")
library("dplyr")
library("limma")
library("matrixStats")
library("cowplot")
library("ggplot2")
library("ggrepel")
library("ggalt")



Project  <- "2018_Turco_Moffett"
Base.dir <- "/Users/rhamilto/Documents/CTR-Manuscripts/2018_Turco_Moffett/2018-Turco-Moffett"
setwd(Base.dir)
list.files(Base.dir)


#source("https://bioconductor.org/biocLite.R")
#biocLite("ggfortify")




message("+-------------------------------------------------------------------------------")
message("+                load in results tables and sample table                        ")
message("+-------------------------------------------------------------------------------")

sample_table          <- read.table("Data_Tables/2018_Turco_Moffett_sample_table.txt", header = T, sep = "\t")
head(sample_table, 30)

norm_matrix           <- read.csv(gzfile("Data_Tables/norm_matrix.csv.gz"), header = T)
rownames(norm_matrix) <- norm_matrix$X
norm_matrix           <- norm_matrix[,-1]
head(norm_matrix)

normCounts_Plac_vs_TOrg <- read.csv(gzfile("Data_Tables/normCounts_Plac_vs_TOrg.csv.gz"), header = T)
normCounts_Plac_vs_TOrg.selection <- normCounts_Plac_vs_TOrg[grepl("(\\bELF1\\b|\\bELF3\\b|\\bELF4\\b|\\bELF5\\b|ERVV-1|GCM1)", normCounts_Plac_vs_TOrg$X),]
normCounts_Plac_vs_TOrg.selection <- normCounts_Plac_vs_TOrg.selection[-c(3),]
colnames(normCounts_Plac_vs_TOrg.selection)[colnames(normCounts_Plac_vs_TOrg.selection)=="TOrg_P043.2"] <- "TOrg_P043_2" 
colnames(normCounts_Plac_vs_TOrg.selection)[colnames(normCounts_Plac_vs_TOrg.selection)=="TOrg_P046.2"] <- "TOrg_P046_2" 
head(normCounts_Plac_vs_TOrg.selection)

normCounts_Str_vs_Plac  <- read.csv(gzfile("Data_Tables/normCounts_Str_vs_Plac.csv.gz"), header = T)
normCounts_Str_vs_Plac.selection  <- normCounts_Str_vs_Plac[grepl("(\\bELF1\\b|\\bELF3\\b|\\bELF4\\b|\\bELF5\\b|ERVV-1|GCM1)", normCounts_Str_vs_Plac$X),]
normCounts_Str_vs_Plac.selection <- normCounts_Str_vs_Plac.selection[-c(3),]
colnames(normCounts_Str_vs_Plac.selection)[colnames(normCounts_Str_vs_Plac.selection)=="TOrg_P043.2"] <- "TOrg_P043_2" 
colnames(normCounts_Str_vs_Plac.selection)[colnames(normCounts_Str_vs_Plac.selection)=="TOrg_P046.2"] <- "TOrg_P046_2"
head(normCounts_Str_vs_Plac.selection)

normCounts_Str_vs_TOrg  <- read.csv(gzfile("Data_Tables/normCounts_Str_vs_TOrg.csv.gz"), header = T)
normCounts_Str_vs_TOrg.selection  <- normCounts_Str_vs_TOrg[grepl("(\\bELF1\\b|\\bELF3\\b|\\bELF4\\b|\\bELF5\\b|ERVV-1|GCM1)", normCounts_Str_vs_TOrg$X),]
normCounts_Str_vs_TOrg.selection <- normCounts_Str_vs_TOrg.selection[-c(3),]
colnames(normCounts_Str_vs_TOrg.selection)[colnames(normCounts_Str_vs_TOrg.selection)=="TOrg_P043.2"] <- "TOrg_P043_2" 
colnames(normCounts_Str_vs_TOrg.selection)[colnames(normCounts_Str_vs_TOrg.selection)=="TOrg_P046.2"] <- "TOrg_P046_2"
head(normCounts_Str_vs_TOrg.selection)

normCounts_Selection           <- as.data.frame(rbindlist(fill=TRUE, use.names=TRUE, 
                                  list(normCounts_Plac_vs_TOrg.selection,normCounts_Str_vs_Plac.selection,
                                       normCounts_Str_vs_TOrg.selection))[,lapply(.SD,mean,na.rm=TRUE), list(X)])
rownames(normCounts_Selection) <- normCounts_Selection$X
normCounts_Selection           <- normCounts_Selection[,-1]

names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "DOrg_M099",   replacement = "DOrg_1")  
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "DOrg_M099",   replacement = "DOrg_1")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "DOrg_M130",   replacement = "DOrg_2")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "DOrg_M134",   replacement = "DOrg_3")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "FOrg_P041",   replacement = "FOrg_1")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "FOrg_P045",   replacement = "FOrg_2")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_085",      replacement = "Pl_1")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_086",      replacement = "Pl_2")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_087",      replacement = "Pl_3")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_241",      replacement = "Pl_4")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_246",      replacement = "Pl_5")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_270",      replacement = "Pl_6")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_278",      replacement = "Pl_7")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Pl_345",      replacement = "Pl_8")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Str_N045",    replacement = "Str_1")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Str_N075",    replacement = "Str_2")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Str_N078",    replacement = "Str_3")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Str_N085",    replacement = "Str_4")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "Str_R031",    replacement = "Str_5")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "TOrg_P041",   replacement = "TOrg_1")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "TOrg_P043_2", replacement = "TOrg_2")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "TOrg_P046_2", replacement = "TOrg_3")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "TOrg_R021",   replacement = "TOrg_4")
names(normCounts_Selection) <- gsub(x = names(normCounts_Selection), pattern = "TOrg_R028",   replacement = "TOrg_5")
head(normCounts_Selection)


message("+-------------------------------------------------------------------------------")
message("+ Customise Heatmaps to give column labels a 45 degree rotation                 ")
message("+-------------------------------------------------------------------------------")
# https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45", ns=asNamespace("pheatmap"))



message("+-------------------------------------------------------------------------------")
message("+ Figure 4E                                                                     ")
message("+-------------------------------------------------------------------------------")

annotation_col = data.frame(Tissue = c("Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi",
                                       "Placental stromal cells", "Placental stromal cells", "Placental stromal cells", "Placental stromal cells", "Placental stromal cells",
                                       "Trophoblast organoids", "Trophoblast organoids", "Trophoblast organoids", "Trophoblast organoids", "Trophoblast organoids"))
rownames(annotation_col) <- colnames(normCounts_Selection)
ann_colors = list(Tissue = c("Placental stromal cells"="springgreen", "Placental villi"="dodgerblue",  "Trophoblast organoids"="lightpink1"))


pdf(paste0(Project, ".Figure.4E.pdf"), onefile=FALSE, width=10, height=5) 
par(bg=NA)
pheatmap(normCounts_Selection,  show_rownames=TRUE, annotation_col = annotation_col, annotation_colors = ann_colors[1],  cutree_row=2, cutree_col=2, cellwidth=25, cellheight=15, fontsize=10, fontsize_row=8) 
dev.off()






message("+-------------------------------------------------------------------------------")
message("+                  get annotations from Ensembl                                 ")
message("+-------------------------------------------------------------------------------")

ensembl                <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host='www.ensembl.org')
#attributes_ensembl <- as.data.frame(listAttributes(ensembl)) 
ensEMBL2id             <- getBM(attributes=c("affy_hugene_2_0_st_v1", "hgnc_symbol", "ensembl_gene_id", "external_gene_name", "entrezgene", "description"),  mart = ensembl) 
ensEMBL2id$description <- gsub("..Source.*", "", ensEMBL2id$description)
ensEMBL2id_probes      <- ensEMBL2id[,c(1,3,4,5)]
ensEMBL2id2            <- ensEMBL2id
ensEMBL2id_probes      <- subset(ensEMBL2id_probes, ensEMBL2id_probes$affy_hugene_2_0_st_v1 != "" ) 

nrow(ensEMBL2id_probes) # 62154 rows
length(unique(ensEMBL2id_probes$ensembl_gene_id)) # 40820 unique gene identifiers
length(unique(ensEMBL2id_probes$affy_hugene_2_0_st_v1)) # 39697 unique affy_hugene_2_0_st_v1 probes 

#interesting.genes.ens <- c("ENSG00000120690", "ENSG00000109381", "ENSG00000163435", "ENSG00000102034", "ENSG00000135374", "ENSG00000137270", "ENSG00000242950")
#interesting.genes.aff <- ensEMBL2id_probes[ensEMBL2id_probes$ensembl_gene_id %in% interesting.genes.ens,]





message("+-------------------------------------------------------------------------------")
message("+                       set up samples and matrix                               ")
message("+-------------------------------------------------------------------------------")

head(norm_matrix)
#colnames(norm_matrix)
#sample_table_original$Sample

colnames(norm_matrix)[colnames(norm_matrix)=="TOrg_P043.2"] <- "TOrg_P043_2" 
colnames(norm_matrix)[colnames(norm_matrix)=="TOrg_P046.2"] <- "TOrg_P046_2" 

#norm_matrix <- norm_matrix[ , order(names(norm_matrix))]

norm_matrix <- norm_matrix[ , c(6:13, 19:21, 23,24, 14:18, 1:3, 22, 4:5)]
head(norm_matrix)

#sample_table_original <-  sample_table[ order(sample_table$Sample_orig_name), ]
#colnames(norm_matrix) == sample_table_original$Sample_orig_name

colnames(norm_matrix) == sample_table$Sample_orig_name
colnames(norm_matrix) <- sample_table$Sample
head(norm_matrix)


rownames(sample_table) <- sample_table$Sample
head(sample_table)

reduced_matrix <- unique(norm_matrix[, colnames(norm_matrix) %in% sample_table$Sample ] )
head(reduced_matrix)

#reduced_matrix <- reduced_matrix[ , order(names(reduced_matrix))]
#sample_table <-  sample_table[ order(sample_table$Sample), ]
#colnames(reduced_matrix) == sample_table$Sample

#sample_table$Tissue <- factor(sample_table$Tissue, levels = c("TOrg", "Str", "Pl"))
#sample_table$Tissue <- factor(sample_table$Tissue, levels = c("TOrg", "Str"))
#Tissue <- sample_table$Tissue




message("+-------------------------------------------------------------------------------")
message("+              filter out probes without gene identifier                        ")
message("+-------------------------------------------------------------------------------")


# filtering probes that have no ensembl gene assigned
probe_stats <- ensEMBL2id_probes %>%
  group_by(affy_hugene_2_0_st_v1) %>%
  summarize(no_of_matches = n_distinct(ensembl_gene_id)) %>%
  filter(no_of_matches > 1)
probe_stats
dim(probe_stats) 
str(probe_stats)
# tells how many probes match more than 1 gene.
# used as a check, but not implemented just yet

# now need to just exclude probes without matching ensembl_gene_id:
probe_to_keep <- ensEMBL2id_probes %>%
  group_by(affy_hugene_2_0_st_v1) %>%
  summarize(no_of_matches = n_distinct(ensembl_gene_id)) %>%
  filter(no_of_matches >= 1)
probe_to_keep
dim(probe_to_keep)

##probe_to_keep <- data.frame(probe_to_keep$affy_hugene_2_0_st_v1)

probe_to_keep <- data.frame(ensEMBL2id_probes$affy_hugene_2_0_st_v1)

probe_to_keep$probe_to_keep <- "probe_to_keep"
colnames(probe_to_keep)[1] <- "affy_hugene_2_0_st_v1"

ensEMBL2id_probes_filtered <- unique(ensEMBL2id_probes[ensEMBL2id_probes$affy_hugene_2_0_st_v1 %in% probe_to_keep$affy_hugene_2_0_st_v1, ])  
ensEMBL2id_filtered <- unique(ensEMBL2id[ensEMBL2id$affy_hugene_2_0_st_v1 %in% probe_to_keep$affy_hugene_2_0_st_v1, ] ) 


# remove probes without matching identifier from counts table as well (call it filtered!):

matrix_filtered <- reduced_matrix[rownames(reduced_matrix) %in% ensEMBL2id_probes_filtered$affy_hugene_2_0_st_v1,  ]  
nrow(matrix_filtered) # 12277
nrow(ensEMBL2id_probes_filtered) # 61643

matrix_filtered_all_samples <- norm_matrix[rownames(norm_matrix) %in% ensEMBL2id_probes_filtered$affy_hugene_2_0_st_v1,  ]  
nrow(matrix_filtered_all_samples) # 12673 
nrow(ensEMBL2id_probes_filtered) # 61643


##matrix_filtered_all_samples <- norm_matrix
##head(matrix_filtered_all_samples)

message("+-------------------------------------------------------------------------------")
message("+                                     QC                                        ")
message("+-------------------------------------------------------------------------------")


#sample_table$Tissue  <- factor(sample_table_original$Tissue, levels = c("TOrg", "Str", "Pl",  "DOrg", "Torg_syn", "FOrg" ))
#sample_table$Colours <- factor(sample_table_original$Colours, levels = c("lightpink1", "springgreen", "dodgerblue",  "mediumpurple", "yellowgreen", "coral" ))

Tissue2 <- sample_table$Tissue
#Tissue2 <- relevel(Tissue2, "Placental_stromal_cells")
design <- model.matrix(~Tissue2, data = matrix_filtered_all_samples)
#design <- model.matrix(~Tissue, data = matrix_filtered)


#levels(Tissue2) <- c("Trophoblast_organoids", "Placental_stromal_cells", "CVS_placenta",  "Decidual_organoids", "Trophoblast_organoids_syncytial", "Fluffy_organoids")


# check clustering with MDS plot:
#pdf(paste(Project, "all_samples", "MDSplot", ".pdf", sep="_"), width=10, height=7, onefile=FALSE)
#par(bg=NA)
plotMDS(matrix_filtered_all_samples, labels=sample_table$samples)
#dev.off()




# PCA plot:
elementTextSize <- 14
#topNum = 500

pca = prcomp(t(matrix_filtered_all_samples))
rv  = rowVars(as.matrix(matrix_filtered_all_samples))
#select = order(rv, decreasing = TRUE)[seq_len(min(topNum, length(rv)))]

pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")



#scores <- data.frame(sample_table, pca$x)

scores <- merge(sample_table, pca$x, by="row.names")


# remove TOrg-P050, FOrg_P045 and FOrg_P041
#scores <- scores[-c(1:5,22),]
scores <- scores[-c(4:5,24),]
scores$Colours <- droplevels(scores$Colours)



plt.4A <- ggplot(scores, aes(x = PC1, y = PC2, fill=factor(scores$Tissue), col=factor(scores$Tissue) )) +
          geom_encircle(alpha = 0.2, show.legend = TRUE) +
          geom_point(size = 5 ) + 
          coord_fixed() +
          xlab(pc1lab) + ylab(pc2lab) + #ggtitle(paste(Project, " PCA Top ", topNum, " MV", sep="")) +
          scale_fill_manual(name="", values=levels(scores$Colours), guide=FALSE) +
          scale_color_manual(name="", labels=c("Trophoblast organoid", "Placental stromal cells", "Placental Villi", "Decidual organoids"), values = levels(scores$Colours)) +
          theme(text=element_text(size=elementTextSize,  family="sans"),
              #  legend.background = element_rect(fill="white", size=.5, color='black', linetype="solid"),   
                legend.position=c(0.05, 0.50), legend.text=element_text(size=elementTextSize), 
                legend.key.height = unit(1.0, "cm")) 
plt.4A


pdf(paste("Turco.Figure.4A.pca.pdf", sep=""), width=6,height=6)
par(bg=NA)
plt.4A
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ Explore PCA Loadings")
message("+-------------------------------------------------------------------------------")

loadings                       <- as.data.frame(pca$rotation)
loadings$affy_hugene_2_0_st_v1 <- rownames(loadings)
loadings                       <- merge(loadings, ensEMBL2id, by="affy_hugene_2_0_st_v1")

topX <- 20

pca.1           <-  loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.topX      <-  pca.1[c(1:topX),]
pca.1.topX.plot <- ggplot(data=pca.1.topX, aes(x=factor(pca.1.topX$external_gene_name,levels=unique(pca.1.topX$external_gene_name)), y=PC1)) + geom_point(size = 3 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.2           <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.topX      <-  pca.2[c(1:topX),]
pca.2.topX.plot <- ggplot(data=pca.2.topX, aes(x=factor(pca.2.topX$external_gene_name,levels=unique(pca.2.topX$external_gene_name)), y=PC2)) + geom_point(size = 3 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



pc.loadings.plt <- plot_grid(pca.1.topX.plot, pca.2.topX.plot, nrow=1, ncol=2 )
pdf(paste("Turco.Figure.4A.loadings.pdf", sep=""), width=10,height=5)
par(bg=NA)
pc.loadings.plt
dev.off()


pdf(paste("Turco.Figure.4A.grid.pdf", sep=""), width=10,height=10)
par(bg=NA)
plot_grid(plt.4A, pc.loadings.plt, nrow=2, ncol=1, rel_heights = c(1, 0.5) )
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                        Prepare av expression table                            ")
message("+-------------------------------------------------------------------------------")


df <- matrix_filtered_all_samples
df <- as.data.frame(t(df))
df$Group <- sample_table_original$Tissue
df2 <- split(df, df$Group)

x <- list()
library(purrr)
for (i in seq_along(df2))  {
  df3 <- data.frame(df2[[i]])
  df3 <- df3[, -ncol(df3)]
  x[[i]] <- map(df3, mean)
  x
}


means_list <- x

levels(df$Group)

for (i in seq_along(levels(df$Group))) {
  names(means_list)[i] <-levels(df$Group)[i]
}

names(means_list)


library(plyr) # also unattach dplyr
# Split List, Apply Function, And Return Results In A Data Frame.
# For each element of a list, apply function then combine results into a data frame.
means_table <- ldply(means_list, function(x){ data.frame(x) })

#rownames(means_table) <- means_table$id
means_tab <- as.data.frame(t(as.data.frame(means_table)))

#means_tab[1,] <- unlist(means_tab[1,])
list1 <- unlist(means_tab[1,])
names(list1) <- list1
colnames(means_tab) <- list1
means_tab <- means_tab[-1,]


#write.csv(means_tab, paste(Project, "means_tab_for_sample_dist_non_tissue_matrix_expr_combat.csv", sep = "_"), quote = F)
#means_tab <- read.csv("CTR_nm390_0001_means_tab_for_sample_dist_non_tissue_matrix_expr_combat.csv", header = T)

write.csv(means_tab, paste(Project, "means_tab_for_sample_dist_matrix_FILTERED_ALL_SAMPLES_expr.csv", sep = "_"), quote = F)
#means_tab <- read.csv("CTR_myt25_0002_means_tab_for_sample_dist_matrix_FILTERED_ALL_SAMPLES_expr.csv", header = T)

#means_tab$X <- gsub("X", "", means_tab$X)
#rownames(means_tab) <- means_tab$X
#means_tab <- means_tab[,-1]





message("+-------------------------------------------------------------------------------")
message("+                          TOrg vs stromal DEGs                                 ")
message("+-------------------------------------------------------------------------------")


#sample_table_original$Tissue <- factor(sample_table_original$Tissue)
#Tissue2 <- sample_table_original$Tissue
Tissue3 <- Tissue2
Tissue3 <- relevel(Tissue3, "Str")
design <- model.matrix(~Tissue3, data = matrix_filtered_all_samples)

#
#
#design <- model.matrix(~Tissue, data = matrix_filtered)
#
#




# DEG analysis with limma



# we compare between different dataset, so add aw into formula:
library(limma)



# fit linear model; From the fit object, calculate empirical Bayes moderated t-statistics.
#fit <- lmFit(matrix_filtered, design) #, weights = aw
fit <- lmFit(matrix_filtered_all_samples, design) #, weights = aw
fit <- eBayes(fit) 

# annotate probes within fit:
fit.df <- as.data.frame(fit)
probes <- data.frame(rownames(fit.df))
colnames(probes)[1] <- "affy_hugene_2_0_st_v1"
probe_anno <- unique(ensEMBL2id)
anno_table <- merge(probes, probe_anno, by.x= T, all.x = T) 
anno_table <- anno_table[order(anno_table$external_gene_name),]
anno_table <- anno_table[!duplicated(anno_table$affy_hugene_2_0_st_v1),]
nrow(anno_table) #  /12277
nrow(fit.df) #  / 12277
fit.df$genes <- anno_table
head(fit.df$genes$external_gene_name)

# check out some comparisons
colnames(design)


fit <- eBayes(fit) 


# get results tables:
topTable(fit, coef="Tissue3TOrg", number = 10)
#topTable(fit, coef="TissueStr", number = 10)
res_TOrg_vs_PlStr <- topTable(fit, coef="Tissue3TOrg", number = Inf)
#res_TOrg_vs_PlStr <- topTable(fit, coef="TissueStr", number = Inf)
sum(res_TOrg_vs_PlStr$adj.P.Val < 0.05 & abs(res_TOrg_vs_PlStr$logFC) > 2)  #     592 (padj 0.05)  /
sum(res_TOrg_vs_PlStr$adj.P.Val < 0.05 & abs(res_TOrg_vs_PlStr$logFC) > 1)  #    2303              / 
sum(res_TOrg_vs_PlStr$adj.P.Val < 0.05 & abs(res_TOrg_vs_PlStr$logFC) > 0)  #    5154              / 

resSig_TOrg_vs_PlStr <- subset(res_TOrg_vs_PlStr, res_TOrg_vs_PlStr$adj.P.Val < 0.05)


resSig_TOrg_vs_PlStr$affy_hugene_2_0_st_v1 <- rownames(resSig_TOrg_vs_PlStr)
resSig_TOrg_vs_PlStr <- merge(resSig_TOrg_vs_PlStr, anno_table, by = "affy_hugene_2_0_st_v1", all.x = T)
resSig_TOrg_vs_PlStr <- resSig_TOrg_vs_PlStr[order(-abs(resSig_TOrg_vs_PlStr$logFC)),]
resSig_TOrg_vs_PlStr <- subset(resSig_TOrg_vs_PlStr, resSig_TOrg_vs_PlStr$adj.P.Val < 0.05)
resSig_TOrg_vs_PlStr_l2fc1 <- subset(resSig_TOrg_vs_PlStr, resSig_TOrg_vs_PlStr$adj.P.Val < 0.05 & abs(resSig_TOrg_vs_PlStr$logFC) > 1)



res_TOrg_vs_PlStr$affy_hugene_2_0_st_v1 <- rownames(res_TOrg_vs_PlStr)
res_TOrg_vs_PlStr <- merge(res_TOrg_vs_PlStr, anno_table, by = "affy_hugene_2_0_st_v1", all.x = T)


write.csv(res_TOrg_vs_PlStr, paste(Project, "res_TOrg_vs_PlStr", "only_TOrg_Str_samples_in_matrix.csv", sep = "_"))
write.csv(resSig_TOrg_vs_PlStr, paste(Project, "resSig_TOrg_vs_PlStr", "only_TOrg_Str_samples_in_matrix.csv", sep = "_"))



#------------------------#------------------------#------------------------#------------------------#------------------------
#      heatmap of TOrg signature genes vs stromal samples ---> all both UP and DOWN now!!!
#------------------------#------------------------#------------------------#------------------------#------------------------

# https://stackoverflow.com/questions/41628450/r-pheatmap-change-annotation-colors-and-prevent-graphics-window-from-popping-up

# for up to 9 colours::::
#mycolors2 <- list(category = brewer.pal(9, "Set1")[1:5])
#names(mycolors2$category) <- levels(annotdf$category)
#cols <- colorRampPalette(brewer.pal(5, "Set1")); mycols <- cols(length(unique(annotdf$category))); mycolors <- list(category=mycols); names(mycolors$category) <- levels(annotdf$category) 


#df <- as.data.frame(sample_table_original[,c("Tissue")])
#colnames(df) <- "Tissue" # "Gestation" # condition
#rownames(df) <- sample_table_original$Sample
#head(df)

logFoldChanceCutOff <- 0.6

#rownames(res.ann.ranked) <- res.ann.ranked$ensembl_gene_id
head(resSig_TOrg_vs_PlStr)
rownames(resSig_TOrg_vs_PlStr) <- resSig_TOrg_vs_PlStr$affy_hugene_2_0_st_v1

nrow(resSig_TOrg_vs_PlStr)

lxl <- rownames( subset(resSig_TOrg_vs_PlStr, abs(resSig_TOrg_vs_PlStr$logFC) > logFoldChanceCutOff))
genes2plot <- unique(lxl)

rows       <- match(genes2plot, row.names(matrix_filtered_all_samples))
mat        <- matrix_filtered_all_samples[rows,]
head(mat)
mat2       <- mat
mat2.df    <- data.frame(affy_hugene_2_0_st_v1=rownames(mat2),mat2)
mat2.ann   <- merge(mat2.df, ensEMBL2id[,-6], by= "affy_hugene_2_0_st_v1" )

mat2.ann$id <- paste(mat2.ann$affy_hugene_2_0_st_v1, mat2.ann$external_gene_name)
head(mat2.ann)

mat2.ann$dups <- duplicated(mat2.ann$id )
mat2.deduped <- subset(mat2.ann, mat2.ann$dups==FALSE)
head(mat2.deduped)
colnames(mat2.deduped)
rownames(mat2.deduped) <- mat2.deduped$id
mat2.new <- mat2.deduped[, -c(1,26:31) ]
head(mat2.new)

# select only stromal vs TORG samples
colnames(mat2.new)
#mat2.new <- mat2.new[, -c(1:5, 6:13, 22)] # to keep only stroma and TOrg
mat2.new <- mat2.new[, -c(1:5, 22)] # to keep stroma, TOrg, and placentaCVS
head(mat2.new)

#annotation_col = data.frame(Tissue = sample_table_original[-c(1:5, 6:13, 22), c(4)])
#annotation_col = data.frame(Tissue = sample_table_original[-c(1:5, 22), c(4)])
annotation_col = data.frame(Tissue = c("Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi", "Placental villi",
                                       "Placental stromal cells", "Placental stromal cells", "Placental stromal cells", "Placental stromal cells", "Placental stromal cells",
                                       "Trophoblast organoids", "Trophoblast organoids", "Trophoblast organoids", "Trophoblast organoids", "Trophoblast organoids"))
rownames(annotation_col) <- colnames(mat2.new)
#ann_colors = list(Tissue = c(Str = "springgreen", TOrg= "lightpink1"))
ann_colors = list(Tissue = c(`Placental stromal cells` = "springgreen", `Trophoblast organoids`= "lightpink1", "Placental villi"="dodgerblue"))

library("pheatmap")

#pdf(paste("Turco.Figure.4C.pdf", sep=""), onefile=FALSE, width=10, height=7) 
#par(bg=NA)
pheatmap(mat2.new,  annotation_col = annotation_col, annotation_colors = ann_colors[1],  cellwidth=25, fontsize=10, fontsize_row=8, show_rownames=F) # , main=paste("Trophoblast genes ::: DEGs log2FoldChange >= ", logFoldChanceCutOff, sep="")
#dev.off()






select.heatmap    <- resSig_TOrg_vs_PlStr
interesting.genes <- c("ENSG00000120690", "ENSG00000109381", "ENSG00000163435", "ENSG00000102034", "ENSG00000135374", "ENSG00000137270", "ENSG00000242950")

select.heatmap[select.heatmap$ensembl_gene_id %in% interesting.genes,]


head(select.heatmap)







# Take the top 50 

rownames(resSig_TOrg_vs_PlStr) <- resSig_TOrg_vs_PlStr$affy_hugene_2_0_st_v1
resSig_TOrg_vs_PlStr2 <- resSig_TOrg_vs_PlStr[order(resSig_TOrg_vs_PlStr$external_gene_name, -abs(resSig_TOrg_vs_PlStr$logFC) ), ]   # sort by id and reverse of abs(value)
resSig_TOrg_vs_PlStr2 <- resSig_TOrg_vs_PlStr2[!duplicated(resSig_TOrg_vs_PlStr2$external_gene_name), ]                       # take the first row within each id
head(resSig_TOrg_vs_PlStr2)
resSig_TOrg_vs_PlStr2 <- subset(resSig_TOrg_vs_PlStr2, resSig_TOrg_vs_PlStr2$adj.P.Val < 0.01)
resSig_TOrg_vs_PlStr2 <- resSig_TOrg_vs_PlStr2[order(abs(resSig_TOrg_vs_PlStr2$logFC), decreasing = T), ]


l30l <- rownames( head(resSig_TOrg_vs_PlStr2[   order( -(resSig_TOrg_vs_PlStr2$logFC)   ), ], n=50) )
#l30l <- rownames( head(resSig_TOrg_vs_PlStr2[   order( -abs(resSig_TOrg_vs_PlStr2$logFC)   ), ], n=50) )
#l30l <- rownames( head(resSig_TOrg_vs_PlStr2[   order( (resSig_TOrg_vs_PlStr2$logFC)   ), ], n=100) )
head(l30l)
length(l30l)
genes2plot <- unique( l30l )
length(genes2plot)




#resSig_TOrg_vs_PlStr2 <- resSig_TOrg_vs_PlStr
#resSig_TOrg_vs_PlStr2 <- resSig_TOrg_vs_PlStr2[order(resSig_TOrg_vs_PlStr2$external_gene_name),]
#resSig_TOrg_vs_PlStr2 <- resSig_TOrg_vs_PlStr2[!duplicated(resSig_TOrg_vs_PlStr2$affy_hugene_2_0_st_v1), ]               
#rownames(resSig_TOrg_vs_PlStr2) <- resSig_TOrg_vs_PlStr2$affy_hugene_2_0_st_v1
#resSig_TOrg_vs_PlStr2 <- subset(resSig_TOrg_vs_PlStr2, resSig_TOrg_vs_PlStr2$adj.P.Val < 0.005)
#l30l <- rownames( head(resSig_TOrg_vs_PlStr2[   order( (resSig_TOrg_vs_PlStr2$logFC)   ), ], n=50) )
#genes2plot <- unique( l30l )





#rows       <- match(genes2plot, row.names(matrix_filtered))
rows       <- match(genes2plot, row.names(matrix_filtered_all_samples))
#mat        <- matrix_filtered[rows,]
mat        <- matrix_filtered_all_samples[rows,]

head(mat)
mat2       <- mat
mat2.df    <- data.frame(affy_hugene_2_0_st_v1=rownames(mat2),mat2)
mat2.ann   <- merge(mat2.df, ensEMBL2id[,-6], by= "affy_hugene_2_0_st_v1" )

mat2.ann$id <- paste(mat2.ann$external_gene_name)
head(mat2.ann)

mat2.ann$dups <- duplicated(mat2.ann$id )
mat2.deduped <- subset(mat2.ann, mat2.ann$dups==FALSE)
head(mat2.deduped)
colnames(mat2.deduped)
rownames(mat2.deduped) <- mat2.deduped$id
mat2.new <- mat2.deduped[, -c(1,26:31) ]
head(mat2.new)

# select only stromal vs TORG samples
colnames(mat2.new)
###mat2.new <- mat2.new[, -c(1:5, 6:13, 22)] # to keep only stroma and TOrg
mat2.new <- mat2.new[, -c(1:5, 22)] # to keep stroma, TOrg, and placentaCVS
head(mat2.new)

#mat2.new <- mat2.deduped[, -c(1, 12:17)]
#annotation_col = data.frame(Tissue = sample_table[, c(4)])
#annotation_col$Tissue <- gsub("Str", "Placental stromal cells", annotation_col$Tissue )
#annotation_col$Tissue <- gsub("TOrg", "Trophoblast organoids", annotation_col$Tissue )
#rownames(annotation_col) <- colnames(mat2.new)
#ann_colors = list(Tissue = c(Str = "springgreen", TOrg= "lightpink1"))


#annotation_col = data.frame(Tissue = sample_table_original[-c(1:5, 6:13, 22), c(4)])
annotation_col = data.frame(Tissue = sample_table_original[-c(1:5, 22), c(4)])
annotation_col$Tissue <- gsub("Pl",   "Placental villi",         annotation_col$Tissue)
annotation_col$Tissue <- gsub("Str",  "Placental stromal cells", annotation_col$Tissue)
annotation_col$Tissue <- gsub("TOrg", "Trophoblast organoids",   annotation_col$Tissue)
head(annotation_col)

rownames(annotation_col) <- colnames(mat2.new)
ann_colors = list(Tissue = c(`Placental stromal cells` = "springgreen", `Trophoblast organoids`= "lightpink1", "Placental villi"="dodgerblue"))
#ann_colors = list(Tissue = c(`Str` = "springgreen", `TOrg`= "lightpink1", `Pl`="dodgerblue"))
head(ann_colors)



pdf("Turco.Figure.4C.pdf", onefile=FALSE, width=12, height=14)
par(bg=NA)
pheatmap(mat2.new, fontsize=13, annotation_col = annotation_col, annotation_colors = ann_colors[1], gaps_col=c(8,5), cellwidth=25, cellheight=15, cutree_rows=5, fontsize_row=13, show_rownames=TRUE)
dev.off()















message("+-------------------------------------------------------------------------------")
message("+-------------------------------------------------------------------------------")
message("+                       pathway analysis/ gene ontology                         ")
message("+-------------------------------------------------------------------------------")
message("+-------------------------------------------------------------------------------")


RESULTS <- resSig_Stroma_vs_TOrg_anno_TORG_SIGNATURE
TITLE <- deparse(substitute(TORG_SIGNATURE_resSig_Stroma_vs_TOrg))
RESULTS$logFC <- abs(RESULTS$logFC)


message("+-------------------------------------------------------------------------------")
message("+                             KEGG pathways                                     ")
message("+-------------------------------------------------------------------------------")

library(pathview)
library(gage)
library(gageData)
library(dplyr)
library(cowplot)

data(korg)
head(korg[,1:3], n=20)
data(kegg.sets.hs)
data(go.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 1)
#data(kegg.sets.mm)
#data(go.sets.mm)
#data(sigmet.idx.mm)
#kegg.sets.mm = kegg.sets.hs[sigmet.idx.mm]
#head(kegg.sets.mm, 1)




Kegg_genes <- RESULTS[!is.na(RESULTS$entrezgene),]
head(Kegg_genes)
nrow(Kegg_genes)
Kegg_genes <- Kegg_genes[order(Kegg_genes$entrezgene, -abs(Kegg_genes$logFC) ), ]        # sort by id and reverse of abs(value)
Kegg_genes <- Kegg_genes[!duplicated(Kegg_genes$entrezgene), ]                       # take the first row within each id
colnames(Kegg_genes)
colnames(Kegg_genes)[colnames(Kegg_genes)=="entrezgene"] <- "Gene_ID" # exchange column name "entrez" to "Gene_ID"

foldchanges = Kegg_genes$logFC
names(foldchanges) = Kegg_genes$Gene_ID
head(foldchanges)

#sum(is.infinite(foldchanges))  # there are some infinite numbers, if use DESeq2, no such problem.
#foldchanges[foldchanges>10]=10
#foldchanges[foldchanges<-10]=-10

colnames(RESULTS)
ResSig <- RESULTS
#ResSig <- subset(ResSig, ResSig$logFC > 1)        # select UP genes!!! (for cell identity)
ResSig <- ResSig[order(ResSig$external_gene_name, -abs(ResSig$logFC) ), ]        # sort by id and reverse of abs(value)
ResSig <- ResSig[!duplicated(ResSig$external_gene_name), ]                       # take the first row within each id

foldchanges_tab <- ResSig[,c("logFC", "entrezgene")]
foldchanges_tab <- na.omit(foldchanges_tab)
foldchanges2 <- foldchanges_tab$logFC
names(foldchanges2) = foldchanges_tab$entrezgene
head(foldchanges2)


keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE, ref=NULL, samp=NULL)
lapply(keggres, head)

data(kegg.gs)
data(go.gs)

keggres.df <- as.data.frame(keggres, row.names=NULL, optional = FALSE, sep= "\t", stringsAsFactors = TRUE)
keggres.df$pathways <- rownames(keggres.df)
keggres.df$pathways <- gsub(" .*", "", keggres.df$pathways)
colnames(keggres.df)
keggres.df = keggres.df[, c(15, 1:14)]
keggres_sig <- sigGeneSet(keggres)
keggsig_gr <- as.data.frame(keggres_sig$greater)
keggsig_less <- as.data.frame(keggres_sig$less)
keggres_sig <- rbind(keggsig_less, keggsig_gr)
write.table(as.data.frame(keggres_sig), file = paste(Project, "Keggres_sig", deparse(substitute(PLACENTAL_SIG_l2fc1_padj0.05)), ".txt", sep = "_" ), sep = "\t")
keggres_sig <- sigGeneSet(keggres)

keggres_sig$greater
keggres_sig$less



# Get the pathways
library(dplyr)
keggrespathways = data.frame(id=rownames(keggres_sig$greater), keggres_sig$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=8) %>% 
  .$id %>% 
  as.character() 
keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids


# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=exp.fc, pathway.id=pid, species="mmu", new.signature=TRUE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

#plot_pathway("hsa04514")








message("+-------------------------------------------------------------------------------")
message("+                                  enrichKegg                                   ")
message("+-------------------------------------------------------------------------------")
library(DOSE)
library(clusterProfiler)
data(geneList)

# organism= "hsa" "mmu" or orther species:  http://www.genome.jp/kegg/catalog/org_list.html
kk <- enrichKEGG(Kegg_genes$Gene_ID, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1) #, readable=TRUE)
head(summary(kk))

str(kk)
kk_results <- as.data.frame(kk)
write.csv(kk_results, file = paste(Project,  TITLE ,"enrichKEGG_", ".csv", sep = "_" ))

# Get the IDs.
keggresids = kk_results$ID
keggresids
#  [1] "hsa03010" "hsa03013" "hsa04151" "hsa05222" "hsa05166" "hsa04217" "hsa04933" "hsa05140" "hsa05164" "hsa05145"


# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=exp.fc, pathway.id=pid, species="human", new.signature=FALSE)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="human"))

plot_pathway_hsa04070 <- pathview(gene.data = foldchanges, pathway.id ="04070", species = "mmu", gene.idtype = "entrez") # kegg.native = FALSE) # out.suffix = "RHMX_04141" , new.signature = FALSE)












message("+-------------------------------------------------------------------------------")
message("+                              Cluster profiler                                 ")
message("+-------------------------------------------------------------------------------")

# predefine gene and genelist:
geneList <- foldchanges2
geneList <- sort(geneList, decreasing = T )
gene <- names(geneList)[abs(geneList) > 2]
gene <- unique(gene)

# GO over-representation test
#ego_bp <- enrichGO(gene = gene, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)
#ego_mf <- enrichGO(gene = gene, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "MF",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)
#ego_cc <- enrichGO(gene = gene, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "CC",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)



# GO over-representation test (gene = all DEGs)
# Gene ID can be mapped to gene Symbol by using paramter readable=TRUE or setReadable function.
ego2_mf <- enrichGO(gene = unique(RESULTS$ensembl_gene_id), OrgDb = org.Hs.eg.db, keytype = 'ENSEMBL', ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = T)
ego2_bp <- enrichGO(gene = unique(RESULTS$ensembl_gene_id), OrgDb = org.Hs.eg.db,keytype = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = T)
ego2_cc <- enrichGO(gene = unique(RESULTS$ensembl_gene_id), OrgDb = org.Hs.eg.db,keytype = 'ENSEMBL', ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = T)
ego2_bp2 <- simplify(ego2_bp, cutoff=0.8, by="p.adjust", select_fun=min)
ego2_cc2 <- simplify(ego2_cc, cutoff=0.8, by="p.adjust", select_fun=min)
ego2_mf2 <- simplify(ego2_mf, cutoff=0.8, by="p.adjust", select_fun=min)
enrichMap(ego2_bp2, n = 50)
enrichMap(ego2_cc2, n = 50)
enrichMap(ego2_mf2, n = 50)

ego2_cc <- as.data.frame(ego2_cc)
ego2_bp <- as.data.frame(ego2_bp)
ego2_mf <- as.data.frame(ego2_mf)

ego2_cc$GO <- "CellularComponent"
ego2_bp$GO <- "BiologicalProcess"
ego2_mf$GO <- "MolecularFuction"

ego2 <- rbind(ego2_bp, ego2_cc, ego2_mf)
head(ego2)
write.csv(ego2, paste(Project, TITLE, "_ego2", "_GO_overrepres_allDEGs_",".csv", sep = ""))





# GO Gene Set Enrichment Analysis
ego3 <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, minGSSize= 10, maxGSSize = 500,  pvalueCutoff = 0.05,verbose  = FALSE)
head(ego3)


# KEGG over-representation test
kk <- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.05)
head(kk)


# KEGG Gene Set Enrichment Analysis
kk2 <- gseKEGG(geneList = geneList, organism = 'hsa', nPerm = 1000, minGSSize = 10, pvalueCutoff = 0.05,verbose = FALSE)
head(kk2)


# KEGG Module over-representation test
mkk <- enrichMKEGG(gene = gene, organism = 'hsa')
as.data.frame(mkk)

# KEGG Module Gene Set Enrichment Analysis
mkk2 <- gseMKEGG(geneList = geneList, minGSSize = 5,  organism = 'hsa', keyType = "kegg")









message("+-------------------------------------------------------------------------------")
message("+                            Gene Ontology (GO)                                 ")
message("+-------------------------------------------------------------------------------")

library(gage)
library(gageData)

data(go.sets.hs)
data(go.subs.hs) 

#Note: "BP", "CC" and "MF", corresponding to biological process, cellular component and molecular function subtrees. It may be more desirable to conduct separated GO enrichment test on each of these 3 subtrees as shown in the example code

gobpsets = go.sets.hs[go.subs.hs$BP] #Biological Process
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)

goccsets = go.sets.hs[go.subs.hs$CC] # cellular component
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)

gomfsets = go.sets.hs[go.subs.hs$MF] # molecular function
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE)
lapply(gomfres, head)







######################     BP - biological process     ########################

gobpres_greater <- as.data.frame(gobpres$greater)
head(gobpres_greater,6)
gobpres_greater$direction <- "Trophoblast_Organoid"   
gobpres_greater$q.val <- as.factor(gobpres_greater$q.val)

gobpres_less <- as.data.frame(gobpres$less)
head(gobpres_less,6)
gobpres_less$direction <- "Placenta_CVS"

gobpres_table <- rbind(subset(gobpres_less, gobpres_less$q.val < 0.1), subset(gobpres_greater, gobpres_greater$q.val < 0.1))
gobpres_less$q.val <- as.factor(gobpres_less$q.val)

gobpres_table$pathways <- rownames(gobpres_table)
colnames(gobpres_table)
gobpres_table = gobpres_table[, c(8, 1:7)]
colnames(gobpres_table)
gobpres_table_top20 <- head( gobpres_table[ order( gobpres_table$q.val ), ], 20)
head(gobpres_table_top20)
gobpres_table_top20$q.val <- as.numeric(gobpres_table_top20$q.val)
gobpres_table_top20$log_q.val <- log10(gobpres_table_top20$q.val)


library(ggplot2) 
#library(plyr) 

# Basic barplot
#gobpres_table_top10_df$pathways <- factor(gobpres_table_top10_df$pathways, levels = gobpres_table_top10_df$pathways[order(gobpres_table_top10_df$greater.q.val)]) 

p <- ggplot(data=gobpres_table_top20, aes(x=reorder(pathways, -q.val), y=-log_q.val, fill=direction)) + 
  geom_bar(stat="identity", alpha = .75) + 
  scale_fill_manual("GO enriched in:", values = c("green4","skyblue2")) + 
  xlab("GO terms") 
p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#+ scale_x_continuous((limits = c(0.0001, 0.05)))  # Horizontal bar plot 
# problem- iether keeing line just before plot to make data,frame and thus continuous values, or without , where graph works but with DISCRETE values.


pdf(paste(Project, "gobpres", "padj0.05", ".pdf", sep="_"), onefile=FALSE, width=10, height=5) 
par(bg=NA)
p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()




######################      MF - moleclar function     ########################



gomfres_greater <- as.data.frame(gomfres$greater)
head(gomfres_greater,6)
gomfres_greater$direction <- "Trophoblast_Organoid"   

gomfres_less <- as.data.frame(gomfres$less)
head(gomfres_less,6)
gomfres_less$direction <- "Placenta_CVS"

gomfres_table <- rbind(subset(gomfres_less, gomfres_less$q.val < 0.03), subset(gomfres_greater, gomfres_greater$q.val < 0.03))
gomfres_less$q.val <- as.factor(gomfres_less$q.val)

gomfres_table$pathways <- rownames(gomfres_table)
colnames(gomfres_table)
gomfres_table = gomfres_table[, c(8, 1:7)]
colnames(gomfres_table)
gomfres_table_top20 <- head( gomfres_table[ order( gomfres_table$q.val ), ], 20)
head(gomfres_table_top20)
gomfres_table_top20$q.val <- as.numeric(gomfres_table_top20$q.val)
gomfres_table_top20$log_q.val <- log10(gomfres_table_top20$q.val)


p<-ggplot(data=gomfres_table_top20, aes(x=reorder(pathways, -q.val), y=-log_q.val, fill=direction)) + 
  geom_bar(stat="identity", alpha = .75) + 
  scale_fill_manual("GO enriched in:", values = c("green4","skyblue2")) + 
  xlab("GO terms") 

p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pdf(paste(Project, "gomfres", "padj0.05", ".pdf", sep="_"), onefile=FALSE, width=10, height=5) 
par(bg=NA)
p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()



######################     CC - cellular component     ########################




goccres_greater <- as.data.frame(goccres$greater)
head(goccres_greater,6)
goccres_greater$direction <- "Trophoblast_Organoid"   

goccres_less <- as.data.frame(goccres$less)
head(goccres_less,6)
goccres_less$direction <- "Placenta_CVS"

goccres_table <- rbind(subset(goccres_less, goccres_less$q.val < 0.05), subset(goccres_greater, goccres_greater$q.val < 0.05))
goccres_less$q.val <- as.factor(goccres_less$q.val)


goccres_table$pathways <- rownames(goccres_table)
colnames(goccres_table)
goccres_table = goccres_table[, c(8, 1:7)]
colnames(goccres_table)
goccres_table_top20 <- head( goccres_table[ order( goccres_table$q.val ), ], 20)
head(goccres_table_top20)
goccres_table_top20$q.val <- as.numeric(goccres_table_top20$q.val)
goccres_table_top20$log_q.val <- log10(goccres_table_top20$q.val)


p<-ggplot(data=goccres_table_top20, aes(x=reorder(pathways, -q.val), y=-log_q.val, fill=direction)) + 
  geom_bar(stat="identity", alpha = .75) + 
  scale_fill_manual("GO enriched in:", values = c("green4","skyblue2")) + 
  xlab("GO terms") 
p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


pdf(paste(Project, "goccres", "padj0.05", ".pdf", sep="_"), onefile=FALSE, width=13, height=5) 
par(bg=NA)
p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()








message("+-------------------------------------------------------------------------------")
message("+                           ReactomePA                                 ")
message("+-------------------------------------------------------------------------------")


logFoldChanceCutOff <- 1 # 
#source("http://bioconductor.org/biocLite.R")
#biocLite('ReactomePA')

#https://bioconductor.org/packages/devel/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
library(ReactomePA)
gene_List <- subset(RESULTS, abs(RESULTS$logFC) > logFoldChanceCutOff)
de <- gene_List$entrezgene
head(de)
x <- enrichPathway(gene=de, pvalueCutoff=0.05, readable=T)
head(as.data.frame(x),20)

png(paste(Project, TITLE, "l2fc1_padj0.05", "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff, "barplot.png", sep="_"), width = 1000, height = 600 )
par(bg=NA)
barplot(x, showCategory=8)
dev.off()

png(paste(Project, TITLE, "l2fc1_padj0.05", "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff,"dotplot.png", sep="_"), width = 1000, height = 600 )
par(bg=NA)
dotplot(x, showCategory=15)
dev.off()

png(paste(Project, TITLE, "l2fc1_padj0.05", "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff,"enrichMap.png", sep="_"), width = 1200, height = 1000 )
par(bg=NA)
cpl <- enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
dev.off()

png(paste(Project, TITLE, "l2fc1_padj0.05", "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff,"cnetplot.png", sep="_"), width = 1500, height = 1500 )
par(bg=NA)
cnetplot(x, categorySize="pvalue", foldChange=geneList)
dev.off()


foldchanges_ordered <- sort(foldchanges, decreasing = T)
foldchanges_ordered <- foldchanges_ordered[!duplicated(names(foldchanges_ordered))]
tab <- as.data.frame(table(names(foldchanges_ordered)))

# Gene Set Enrichment Analysis
y <- gsePathway(foldchanges_ordered, nPerm=1000,
                minGSSize=5, pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE)
ress <- as.data.frame(y)
head(ress)
pathways <- ress$Description

write.csv(as.data.frame(y), paste(Project,  TITLE, "l2fc1_padj0.05", "gsePathways_reactomeRA","logFoldChanceCutOff", logFoldChanceCutOff, ".csv", sep = ""))

png(paste(Project, TITLE, "l2fc1_padj0.05", "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff, "GSEAenrichMap.png", sep="_"), width = 1500, height = 1500 )
par(bg=NA)
enrichMap(y)
dev.off()












#------------------------#------------------------#------------------------#------------------------#------------------------
#      heatmap of TOrg signature genes vs stromal samples
#------------------------#------------------------#------------------------#------------------------#------------------------


annotation_col = data.frame(Tissue = sample_table_original[-c(1:5, 22), c(4)])
rownames(annotation_col) <- colnames(mat2.new)
ann_colors = list(Tissue = c(Str = "springgreen", TOrg= "lightpink1", Pl= "dodgerblue"))
head(annotation_col)
head(ann_colors)


logFoldChanceCutOff <- 0.6


#reduced_matrix <- unique(norm_matrix[, colnames(norm_matrix) %in% sample_table$Sample ] )

lxl <- unique(ensEMBL2id[ ensEMBL2id$external_gene_name %in% ESC_genes$external_gene_name , ] )
lxl <- lxl[!is.na(lxl$affy_hugene_2_0_st_v1),]
genes2plot <- unique(lxl$affy_hugene_2_0_st_v1)

lxl <- lxl[order(lxl$external_gene_name),]


rows       <- match(genes2plot, row.names(matrix_filtered_all_samples))
mat        <- matrix_filtered_all_samples[rows,]
mat        <- na.omit(mat)
head(mat, 10)
mat2       <- mat
mat2.df    <- data.frame(affy_hugene_2_0_st_v1=rownames(mat2),mat2)
mat2.ann   <- merge(mat2.df, ensEMBL2id[,-6], by= "affy_hugene_2_0_st_v1" )

mat2.ann$id <- paste(mat2.ann$affy_hugene_2_0_st_v1, mat2.ann$external_gene_name)
mat2.ann$dups <- duplicated(mat2.ann$id )
mat2.deduped <- subset(mat2.ann, mat2.ann$dups==FALSE)
head(mat2.deduped, 20)
colnames(mat2.deduped)

mat2.deduped$meanTORG <- rowMeans(subset(mat2.deduped, select = c(20:22,24:25)), na.rm = TRUE)
#mat2.deduped$meanPL <- rowMeans(subset(mat2.deduped, select = c(7:14)), na.rm = TRUE)

mat2.deduped <- mat2.deduped[order(mat2.deduped$meanTORG, decreasing = T),]
#mat2.deduped <- mat2.deduped[order(mat2.deduped$meanPL, decreasing = T),]

#mat2.deduped$dups <- paste(mat2.deduped$meanTORG, mat2.deduped$external_gene_name)
mat2.deduped$dups <- paste(mat2.deduped$external_gene_name)

mat2.deduped <- mat2.deduped[!duplicated(mat2.deduped$dups),]
rownames(mat2.deduped) <- mat2.deduped$external_gene_name


mat3 <- mat2.deduped[, c(27:30,32)] # TORG only
#mat4 <- mat2.deduped[, c(7:14)]

#rownames(mat2.deduped) <- mat2.deduped$id
mat2.new <- mat2.deduped[, -c(1,26:32) ]
head(mat2.new)

# select only stromal vs TORG samples
colnames(mat2.new)
#mat2.new <- mat2.new[, -c(1:5, 6:13, 22)] # to keep only stroma and TOrg
mat2.new <- mat2.new[, -c(1:5, 22)] # to keep stroma, TOrg, and placentaCVS
mat2.new <- mat2.new[, -c(1:5, 6:13, 14:18, 22)] # to keep TOrg
head(mat2.new)



pdf(paste(Project, "_mouse_TSC-enriched_genes__Stroma_PlacCVS_TOrg_CountMatrixHeatmap", "_clust_by_Pl.pdf", sep=""), onefile=FALSE, width=7, height=7) 
par(bg=NA)
pheatmap(mat2.new, fontsize=10, annotation_col = annotation_col, annotation_colors = ann_colors[1], fontsize_row=8, show_rownames=F, cluster_rows = F, cluster_cols = TRUE) 
dev.off()

pdf(paste(Project, "_mouse_TSC-enriched_genes_Stroma_PlacCVS_TOrg_CountMatrixHeatmap_top50", "_clust_by_Pl.pdf", sep=""), onefile=FALSE, width=7, height=10) 
par(bg=NA)
pheatmap(mat2.new[c(1:50),], fontsize=10, annotation_col = annotation_col, annotation_colors = ann_colors[1], fontsize_row=8, show_rownames=T, cluster_rows = F, cluster_cols = TRUE) 
dev.off()

pdf(paste(Project, "_mouse_TSC-enriched_genes_Stroma_PlacCVS_TOrg_CountMatrixHeatmap_top100", "_clust_by_Pl.pdf", sep=""), onefile=FALSE, width=8, height=15) 
par(bg=NA)
pheatmap(mat2.new[c(1:100),], fontsize=10, annotation_col = annotation_col, annotation_colors = ann_colors[1], fontsize_row=8, show_rownames=T, cluster_rows = F, cluster_cols = TRUE) 
dev.off()

pdf(paste(Project, "_mouse_TSC-enriched_genes_Stroma_PlacCVS_TOrg_CountMatrixHeatmap_top30", "_clust_by_Pl.pdf", sep=""), onefile=FALSE, width=7, height=7) 
par(bg=NA)
pheatmap(mat2.new[c(1:30),], fontsize=10, annotation_col = annotation_col, annotation_colors = ann_colors[1], fontsize_row=8, show_rownames=T, cluster_rows = F, cluster_cols = TRUE) 
dev.off()







mat3$meanTORG <- as.numeric(mat3$meanTORG)
mat3$rank <- NA
mat3 <- mat3[order( mat3$external_gene_name),] 
mat3$rank[order(mat3$meanTORG,  decreasing = T)] <- 1:nrow(mat3)
mat3$rank <- as.numeric(mat3$rank)
head(mat3)

coef(lm(meanTORG ~ rank, data = mat3))



pdf(paste(Project, "_TSC-enriched genes_rank_expression_TOrg_Scatterplot", ".pdf", sep=""), onefile=FALSE, width=10, height=7) 
par(bg=NA)
ggplot(mat3, aes(x=rank, y=meanTORG)) + 
  geom_point(colour = "black", size = 2, alpha = 0.3) +
  geom_abline(intercept = 7.80394470, slope = -0.01569908, colour = "grey", alpha = 0.5) +
  geom_vline(xintercept = 30, colour = "red", alpha = 0.5) +
  xlab("Rank") +
  ylab("Mean expression in Trophoblast Organoids") +
  ggtitle("Expression distribution of mouse TSC-enriched genes in Trophoblast Organoids")
dev.off()






head(mat2.new)



pdf(paste(Project, "_all_human_like_ESC_genes_Stroma_PlacCVS_TOrg_CountMatrixHeatmap","_top50" , ".pdf", sep=""), onefile=FALSE, width=10, height=10) 
par(bg=NA)
pheatmap(mat2.new[c(1:50),], fontsize=10, annotation_col=df, fontsize_row=8, show_rownames=T, cluster_rows = F, cluster_cols = TRUE) 
dev.off()




# PCA plot:
library(matrixStats)
elementTextSize <- 8
topNum = 50

pca = prcomp(t(mat2.new))
rv = rowVars(as.matrix(mat2.new))
select = order(rv, decreasing = TRUE)[seq_len(min(topNum, length(rv)))]

pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

scores <- data.frame(sample_table_original[-c(1:5, 22),], pca$x, Tissue2[-c(1:5, 22)])



library(cowplot)
library(ggplot2)
library(ggrepel)
pdf(paste(Project, "_TOrg_Str_PlCSV_", "top50_ESC_markers_","PCA", ".pdf", sep=""), width=8,height=6)
par(bg=NA)
ggplot(scores, aes(x = PC1, y = PC2, col = factor(scores$Tissue2) )) +
  geom_point(size = 5 ) + #geom_text_repel(aes(label=scores$Sample), col = "black") +
  xlab(pc1lab) + ylab(pc2lab) + #ggtitle(paste(Project, " PCA Top ", topNum, " MV", sep="")) +
  scale_colour_manual(name="Tissue", values = c( "blue2",  "purple", "green4", "#5f022f", "violetred", "pink", "cornflowerblue", "blue")) +
  theme(text = element_text(size=elementTextSize)) 
dev.off()










message("+-------------------------------------------------------------------------------")
message("+                            Visualize GO terms                                 ")
message("+-------------------------------------------------------------------------------")



library(GOplot)


# Name	          Description	                                                                             Dimension
# 
# EC$eset	        Data frame of normalized expression values                                               20644 x 7
# EC$genelist	    Data frame of differentially expressed genes (adjusted p-value < 0.05)	                 2039 x 7
# EC$david	      Data frame of results from a functional analysis of the differentially expressed genes 	 174 x 5
# EC$genes	      Data frame of selected genes with logFC	                                                 37 x 2
# EC$process	    Character vector of selected enriched biological processes	                             e.g 7


#data(EC)
#str(EC) # list of data frames
#circ <- circle_dat(EC$david, EC$genelist)

which(grepl("collagen metabolic process", ego3$Description)) # 18
ego3$core_enrichment[34]

# Generate the plotting object: 

rownames(RESULTS) <- RESULTS$affy_hugene_2_0_st_v1
colnames(RESULTS)
ResSig <- merge(RESULTS, av_expressions_filtr, by = "row.names" )
ResSig <- ResSig[order(ResSig$external_gene_name, -abs(ResSig$logFC) ), ]        # sort by id and reverse of abs(value)
ResSig <- ResSig[!duplicated(ResSig$external_gene_name), ]                       # take the first row within each id
colnames(ResSig)
ResSig <- unique(ResSig)
genes <- ResSig[,c(11,3)]
genes <- genes[order(genes$external_gene_name, -abs(genes$logFC) ), ]        # sort by id and reverse of abs(value)
genes <- genes[!duplicated(genes$external_gene_name), ]                       # take the first row within each id
colnames(genes) <- c("ID", "logFC")
ResSig2 <- ResSig[, c(11,3,4,5,6,7,8,10,12,14,15)]
colnames(ResSig2) <- c( "ID", "logFC", "AveExpr", "t",  "P.Value"  ,  "adj.P.Val" ,  "B" ,  "ensembl_gene_id",   "entrezgene",  "av_Pl", "av_TOrg" )
ResSig2 <- unique(ResSig2)
ResSig2 <- ResSig2[order(ResSig2$ID, -abs(ResSig2$logFC) ), ]        # sort by id and reverse of abs(value)
ResSig2 <- ResSig2[!duplicated(ResSig2$ID), ]                       # take the first row within each id


head(ego2) #  here ego2 are all sigDEGs with --- padj 0.05 & l2fc1 thresholds
colnames(ego2)[colnames(ego2)=="p.adjust"] <- "adj_pval" # 
colnames(ego2)[colnames(ego2)=="geneID"] <- "genes" # 
colnames(ego2)[colnames(ego2)=="core_enrichment"] <- "genes" # 
colnames(ego2)[colnames(ego2)=="Description"] <- "term" # 
ego2$genes <- gsub("/", ", ", ego2$genes)
colnames(ego2)[colnames(ego2)=="GO"] <- "Category" # 


ResSig2$logFC <- as.numeric(ResSig2$logFC)


process <- c("embryonic placenta development", "gland development")


EC2 <- list(eset = reduced_matrix[,-c(1:8)], genelist = ResSig2, ego2 = ego2, genes = genes, process = process, ego2 = ego2)



head(EC2$genelist, 10)
head(EC2$ego2)
head(EC2$genes)
head(EC2$process)


circ <- circle_dat(EC2$ego2, EC2$genes)
colnames(circ)

circ$zscore
circ$count

#circ2 <- circ[!is.na(circ$zscore),]


# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.9) # circ made from ego2
reduced_circ$logFC <- as.numeric(reduced_circ$logFC)
# now some clean-up:
reduced_circ <- reduced_circ[,]




# ...and plot it
GOBubble(circ, labels = 2)




pdf(paste(Project, "PlacCVS_vs_TOrg", "reduced_circ","GOBar", ".pdf", sep="_"), onefile=FALSE, width=12, height=5) 
par(bg=NA)
GOBar(reduced_circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'))
#GOBar(circ2, display = 'multiple')
dev.off()



# Generate the bubble plot with a label threshold of 2

pdf(paste(Project, TITLE,"reduced_circ","GOBubble", ".pdf", sep="_"), onefile=FALSE, width=20, height=12) 
par(bg=NA)
GOBubble(reduced_circ, labels = 2)
dev.off()

pdf(paste(Project, TITLE,"circ","GOBubble", ".pdf", sep="_"), onefile=FALSE, width=20, height=12) 
par(bg=NA)
GOBubble(circ, labels = 2)
dev.off()







# Generate a circular visualization of the results of gene- annotation enrichment analysis

IDs <- c('GO:0001892', 'GO:0048732')

pdf(paste(Project, TITLE,  "padj0.05", "GOCircle", ".pdf", sep="_"), onefile=FALSE, width=10, height=10) 
par(bg=NA)
GOCircle(circ, nsub = IDs, label.size = 2, zsc.col = c('red', 'white', 'blue'))

dev.off()






#
#
#   junction

process <- c("cell junction organization", "cell-cell junction organization", "apical junction complex", "cell-cell adherens junction", "occluding junction", "ruffle membrane", "basolateral plasma membrane", "cadherin binding", "extrinsic component of membrane", "regulation of small GTPase mediated signal transduction", "apical part of cell" )
chord <- chord_dat(circ, subset(EC2$genes, abs(EC2$genes$logFC) > 1), process)
head(chord)

pdf(paste(Project, TITLE,  "chord_l2fc1", "JUNCTIONS","GOChord", ".pdf", sep="_"), onefile=FALSE, width=10, height=10) 
#pdf(paste(Project, TITLE,  "chord_l2fc1", "JUNCTIONS","GOChord_legend", ".pdf", sep="_"), onefile=FALSE, width=20, height=10) 
par(bg=NA)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4, process.label =15, lfc.min=-2, lfc.max = 2 ,  limit = c(3,4))
dev.off()


#
#
#   sterols

process <- c("alcohol biosynthetic process", "cholesterol biosynthetic process", "secondary alcohol biosynthetic process", "alcohol metabolic process", "sterol biosynthetic process", "glycogen catabolic process", "extrinsic component of membrane", "macromolecule deacylation" )
chord <- chord_dat(circ, subset(EC2$genes, abs(EC2$genes$logFC) > 1), process)
head(chord)

pdf(paste(Project, TITLE,  "chord_l2fc1", "STEROLS","GOChord", ".pdf", sep="_"), onefile=FALSE, width=7, height=7) 
#pdf(paste(Project, TITLE,  "chord_l2fc1", "STEROLS","GOChord_legend", ".pdf", sep="_"), onefile=FALSE, width=20, height=10) 
par(bg=NA)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4, process.label =15, lfc.min=-2, lfc.max = 2 ,  limit = c(3,4))
dev.off()




# testing which terms to go for:
process <- unique(circ$term)
chord <- chord_dat(circ, subset(EC2$genes, abs(EC2$genes$logFC) > 1), process)
head(chord)

pdf(paste(Project, TITLE,  "chord_l2fc1", "all_terms","GOChord_legend", ".pdf", sep="_"), onefile=FALSE, width=35, height=20) 
par(bg=NA)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 7, process.label =15, lfc.min=-2, lfc.max = 2 ,  limit = c(5,7))
dev.off()


pdf(paste(Project,TITLE,  "chord_l2fc1", "all_terms","GOChord", ".pdf", sep="_"), onefile=FALSE, width=15, height=15) 
par(bg=NA)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 7, process.label =15, lfc.min=-2, lfc.max = 2 ,  limit = c(5,7))
dev.off()
















sessionInfo() 


























