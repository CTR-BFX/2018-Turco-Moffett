# Trophoblast organoids model early human placental development #

** Margherita Y. Turco<sup>1,2,3,§</sup>, Lucy Gardner<sup>1,3</sup>, Richard Kay<sup>4</sup>, Malwina Prater<sup>3</sup>, Alasdair McWhinnie<sup>5</sup>, M. Hollinshead<sup>1</sup>, Laura Esposito<sup>1</sup>, Frank Reimann<sup>4</sup>, Fiona Gribble<sup>4</sup>, Andrew Sharkey<sup>1,3</sup>, Steven G.E. Marsh<sup>5</sup>, Stephen O’Rahilly<sup>4</sup>, Myriam Hemberger<sup>3,6</sup>, Graham J. Burton<sup>2,3,7,§</sup> and Ashley Moffett<sup>1,3,7,§</sup> **

<sup>1</sup> Department of Pathology, University of Cambridge, UK
<sup>2</sup> Department of Physiology, Neurobiology and Development, University of Cambridge, UK
<sup>3</sup> Centre for Trophoblast Research, University of Cambridge, UK
<sup>4</sup> Metabolic Research Laboratories, Wellcome Trust-MRC Institute of Metabolic Science, Addenbrooke’s Hospital, Cambridge, UK
<sup>5</sup> Anthony Nolan Research Institute, Royal Free Hospital, London, UK
<sup>6</sup> Epigenetics Programme, The Babraham Institute, Cambridge, UK
<sup>7</sup> Co-last authors
<sup>§</sup> Correspondence: M.Y. Turco (myt25@cam.ac.uk), G.J. Burton (gjb2@cam.ac.uk) and A. Moffett (am486@cam.ac.uk)



### Citation ###

Turco, M.Y., Gardner, L., Kay, R., Prater, M., McWhinnie, A., Hollinshead., M., Esposito, L., Reimann, F., Gribble, F., Sharkey, A., Marsh, S.G.E., O’Rahilly, S., Hemberger, M., Burton, G.J. and Moffett, A. (2018) Trophoblast organoids model early human placental development.

### Abstract ###

To appear on publication

### Data Processing:

The dataset consists of groups: CVS samples (x8), Trophoblast Organoids (x5), Placental stromal cells (x5).

The normalised matrix was prepared using lumi. Microarray probes without gene identifier (ensembl gene id) were filtered out. Initial QC included PCA, MDS plot, etc.   Comparisons were performed using the limma package (3.34.8) with results corrected for multiple testing using False Discovery Rate (FDR) testing. Finally the quality of the data was assessed and the correlation of the samples in the groups compared. Heatmaps were generated using the pheatmap function of the R package 'pheatmap' (1.0.8), which uses euclidean method to obtain the distance matrix and complete agglomeration method for clustering.  For the gene heatmaps, the input is the normalized intensity matrix. GO term enrichment was obtained with R package 'clusterProfiler' (3.6.0) with function enrichGO and chord plots were generated using the R package 'GOplot' (1.0.2).

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


### Sample Table ###

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

### Links ###

Description   | URL
------------- | ----------
Publication   | [Journal](http://) and [DOI](http://) <br>(<i>To be updated on publication</i>)
Raw Data      | ArrayExpress EMBL-EBI [E-MTAB-6683](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6683) <br>(<i>Data to be released on publication</i>)

### Contact

Contact Malwina Prater (mn367@cam.ac.uk) or Russell Hamilton (rsh46@cam.ac.uk) for bioinformatics related queries.
