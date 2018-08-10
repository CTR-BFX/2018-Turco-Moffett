# Long-tern, functional human Trophoblast organoids provide a model for maternal-fetal interactions in early pregnancy #

**Margherita Y. Turco<sup>1,2,3,§</sup>, Lucy Gardner<sup>1,3</sup>, Richard Kay<sup>4</sup>, Russell S. Hamilton<sup>2,3</sup>, Malwina Prater<sup>2,3</sup>, Michael Hollinshead<sup>1</sup>, Alasdair McWhinnie<sup>5</sup>, Laura Esposito<sup>1</sup>, Ridma Fernando<sup>2,3</sup>, Helen Skelton<sup>1</sup>, Frank Reimann<sup>4</sup>, Fiona Gribble<sup>4</sup>, Andrew Sharkey<sup>1,3</sup>, Steven G.E. Marsh<sup>5,6</sup>, Stephen O’Rahilly<sup>4</sup>, Myriam Hemberger<sup>3,7</sup>, Graham J. Burton<sup>2,3,8</sup> and Ashley Moffett<sup>1,3,8,§</sup>**

<sup>1</sup> Department of Pathology, University of Cambridge, UK
<sup>2</sup> Department of Physiology, Neurobiology and Development, University of Cambridge, UK
<sup>3</sup> Centre for Trophoblast Research, University of Cambridge, UK
<sup>4</sup> Metabolic Research Laboratories, Wellcome Trust-MRC Institute of Metabolic Science, Addenbrooke’s Hospital, Cambridge, UK
<sup>5</sup> Anthony Nolan Research Institute, Royal Free Hospital, London, UK
<sup>6</sup>UCL Cancer Institute, Royal Free Campus, London, UK
<sup>7</sup> Epigenetics Programme, The Babraham Institute, Cambridge, UK
<sup>8</sup> Co-last authors
<sup>§</sup> Correspondence: M.Y. Turco (myt25@cam.ac.uk) and A. Moffett (am486@cam.ac.uk)



### Citation ###

Turco, M.Y., Gardner, L., Kay, R., Hamilton, R.S., Prater, M., McWhinnie, A., Hollinshead., M., Esposito, L., Fernando, R., Skelton, H., Reimann, F., Gribble, F., Sharkey, A., Marsh, S.G.E., O’Rahilly, S., Hemberger, M., Burton, G.J. and Moffett, A. (2018) Long-tern, functional human Trophoblast organoids provide a model for maternal-fetal interactions in early pregnancy.

### Abstract ###

To appear on publication

### Microarray Data Processing:

The dataset consists of groups: CVS samples (x8), Trophoblast Organoids (x5), Placental stromal cells (x5).

The normalised matrix was prepared using lumi. Microarray probes without gene identifiers (ensembl gene id) were filtered out. Initial QC included PCA, MDS plot, etc. Principal component loadings were extracted to reveal the top 20 genes contributing to each PC. Comparisons were performed using the limma package (3.34.8) with results corrected for multiple testing using False Discovery Rate (FDR) testing. Finally the quality of the data was assessed and the correlation of the samples in the groups compared. Heatmaps were generated using the pheatmap function of the R package 'pheatmap' (1.0.8), which uses euclidean method to obtain the distance matrix and complete agglomeration method for clustering.  For the gene heatmaps, the input is the normalized intensity matrix. GO term enrichment was obtained with R package 'clusterProfiler' (3.6.0) with function enrichGO and chord plots were generated using the R package 'GOplot' (1.0.2).

Placental identity was derived from comparison Placental stromal cells vs CVS placenta, and the upregulated genes defined placental trophoblast identity. Stromal genes were identified as significantly upregulated (l2fc > 1) compared to Placenta and TOrg using limma linear model and calculating empirical Bayes moderated t-statistics. The threshold of l2fc > 1 was set up low as Placenta samples already contain up to 40% stromal cells which would affect global gene expression levels. For the comparison of TOrg vs Placenta samples (difference between natural and organoid system), stromal genes were removed from the filtered matrix, and limma linear model and calculating empirical Bayes moderated t-statistics were applied. TOrg signatures (from comparison of of TOrg vs stroma) were characterised by pathway/GO analyses. The comparison between TOrg and Placenta samples with mouse ESC signature genes was done using gene list kindly provided by Myriam Hemberger. Pathway analysis was performed using databases: GO, Kegg, and Reactome,  with R packages: dose (3.2.0),  clusterProfiler (3.4.4), gage (2.26.1) and ReactomePA (1.20.2).

Resource       | URL
-------------- | --------------
GRCh38         | [Link](http://mar2016.archive.ensembl.org/index.html)


### Script to reproduce figure 4 ###

Figure    | Output Filename                             | Description  
--------- | ------------------------------------------- | ------------------------
4A        | Turco_Figure.4A.pdf         | PCA and PC loadings plots
4B        | Turco_Figure.4B.pdf         | Heatmap
4C        | Turco_Figure.4C.pdf         | Heatmap (top 50 by l2fc)
4D        | Turco_Figure.4D.pdf         | Chord GO Plot
4E        | Turco_Figure.4E.pdf         | Custom heatmap for selected genes


### Microarray Sample Table ###

Sample_orig_name	| Sample	| Tissue_original	|
------------|--------|------------------------
Pl_085	    | Pl_1	 | Placental villi  
Pl_086	    | Pl_2	 | Placental villi	 
Pl_087	    | Pl_3	 | Placental villi	 
Pl_241	    | Pl_4	 | Placental villi	 
Pl_246	    | Pl_5	 | Placental villi	 
Pl_270	    | Pl_6	 | Placental villi	 
Pl_278    	| Pl_7	 | Placental villi	 
Pl_345	    | Pl_8	 | Placental villi	 
TOrg_P041	  | TOrg_1 | Trophoblast organoids	 
TOrg_P043_2	| TOrg_2 | Trophoblast organoids	 
TOrg_P046_2	| TOrg_3 | Trophoblast organoids	 
TOrg_R021	  | TOrg_4 | Trophoblast organoids	 
TOrg_R028	  | TOrg_5 | Trophoblast organoids	 
Str_N045	  | Str_1	 | Placental stromal cells	 
Str_N075	  | Str_2	 | Placental stromal cells	 
Str_N078	  | Str_3	 | Placental stromal cells	 
Str_N085	  | Str_4	 | Placental stromal cells	 
Str_R031	  | Str_5	 | Placental stromal cells
DOrg_M099   | DOrg_1 | Decidual organoids    
DOrg_M130   | DOrg_2 | Decidual organoids    
DOrg_M134   | DOrg_3 | Decidual organoids    	 

### EPIC Methylation Array Data Processing:
Genomic DNA bisulfite (BS) and oxidative bisulfite (oxBS) conversion were performed using the CEGX TrueMethyl kit (Cambridge Epigenetix / NuGEN ) and used for microarray-based DNA methylation analysis, performed at GenomeScan (GenomeScan B.V., Leiden, The Netherlands) on the HumanMethylation850 BeadChip (Illumina, Inc., San Diego, CA, U.S.A). This array interrogates over 850,000 CpG sites representing about 99% of the RefSeq genes.  The resulting iDAT files were imported and analysed using ChAMP (v2.9.10)[1,2]. Samples were processed filtering for a probe detection p-value <= 0.01, probes with a beadcount <3 in at least 5% of samples, no CpG and known SNPs[3] at probe starts, probes aligning to multiple locations,  and QC using the on array control probes. Of the total probes on the array 755577 passed the filtering and QC steps. The BMIQ[4] method was used to normalisation the two probe types present on the array. Beta methylation values from the EPIC array range from 0 (unmethylated) to 1 (methylated) and are equivalent of percentage methylation. Genomic annotations were imported from FDb.InfiniumMethylation.hg19 and IlluminaHumanMethylationEPICmanifest[5]. LINE1 elements were downloaded as tables from UCSC Genome browser for hg19[6].

Rscript to recreate Methylation figures []()

### EPIC Methylation Array Sample Table ###

Sample_orig_name | Sample  | Tissue_original |
-----------------|---------|-----------------|
Pl_246	         | Pl_5	   | Placental villi	 
Pl_270	         | Pl_6	   | Placental villi	 
Pl_278           | Pl_7	   | Placental villi	 
Pl_345	         | Pl_8	   | Placental villi	 
TOrg_R21         | TOrg_4  | Trophoblast organoids
TOrg_R51         | TOrg_10 | Trophoblast organoids
TOrg_R70         | TOrg_12 | Trophoblast organoids
TOrg_R72         | TOrg_14 | Trophoblast organoids

### References ###

1.	Morris, T. J. et al. ChAMP: 450k Chip Analysis Methylation Pipeline. Bioinformatics 30, 428–430 (2014).
2.	Aryee, M. J. et al. Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics 30, 1363–1369 (2014).
3.	Zhou, W., Laird, P. W. & Shen, H. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Res. 45, e22 (2017).
4.	Teschendorff, A. E. et al. A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data. Bioinformatics 29, 189–196 (2013).
5.	Fortin, J., Triche, T. J. & Hansen, K. D. Preprocessing , normalization and integration of the Illumina HumanMethylationEPIC array.
6.	Karolchik, D. et al. The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 32, D493-6 (2004).



### Links ###

Description   | URL
------------- | ----------
Publication   | [Journal](http://) and [DOI](http://) <br>(<i>To be updated on publication</i>)
Raw Data (microarray)     | ArrayExpress EMBL-EBI [E-MTAB-6683](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6683) <br>(<i>Data to be released on publication</i>)
Raw Data (EPIC)     | [EPIC_IDATs files](EPIC_IDATs) available in this repository <br>ArrayExpress EMBL-EBI [E-MTAB-XXXX](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-XXXX) <br>(<i>Data to be released on publication</i>)
### Contact

Contact Malwina Prater (mn367@cam.ac.uk) or Russell Hamilton (rsh46@cam.ac.uk) for bioinformatics related queries.
