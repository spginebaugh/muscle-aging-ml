## Multicellular Factor Analysis for the Identification of the Transcriptomic Signal of Aging in Muscle Tissue

### Abstract


A common method to identify potential therapeutic targets in transcriptomic data is to use differential expression (DE) to identify highly differentially expressed genes between two conditions. When it comes to single-cell/single-nucleus datasets, DE is great at identifying transcriptomic differences between cell types (a method commonly used to identify marker genes). ScRNAseq is also very useful for identifying differences in cell type quantities between samples/conditions. However, DE is often not ideal for identifying within-cluster, between-sample/condition transcriptomic changes in scRNAseq datasets (e.g., what are the transcriptomic differences between muscle endothelial cells from young and old individuals). This is because batch effects, noise, variation in the number of cells per sample, and variation in the number of samples/donors per condition are often quite large in scRNAseq, leading to in inaccurate results from DE.

Cells can have both cell-type specific and shared responses to the same cue, and well as feedback loops that occur through cell-cell interactions and between-organ signals (e.g. hormones). These multicellular changes are sometimes refered to as "multicellular programs". Between-sample DE looks at one cell-type at a time and does not consider these multicellular programs, which can result in "missing the forest for the trees". This is especially true for aging, which clearly impacts every celltype in the body in numerous ways.

We need more sophisticated methods that are able to extract important signals from noise and identify both celltype-specific and multicellular program changes. One such approach is to use Blind Source Separation techniques, such as PCA, ICA, and NMF. However, when applied directly to scRNAseq data, these methods will mainly identify differences between cell types. This limitation can be overcome by taking our genes x cells matrix, and turn it into a pseudobulk celltype x genes x samples tensor. We can then decompose the tensor with blind source separation to pull out the important biological signals.

Here, we analyze single-nucleus RNAseq data from [Perez et al. 2022](https://doi.org/10.18632/aging.204435) using [MOFAcell](https://doi.org/10.7554/eLife.93161) -- a factor analysis approach to single cell data. Ideally, this method will identify a factor that is highly correlated with age. That "Aging Factor" should then contain a scored list of how much each gene contributes to aging in each cell type, independent of technical noise or other biological factors.

***

### Analyses & Scripts

The code to for the preprocessing of the dataset can be found [here](https://github.com/spginebaugh/muscle_aging_ML/tree/main/scripts/preparation)

The results and visulazation of the analysis can be found [here](https://github.com/spginebaugh/muscle_aging_ML/blob/main/scripts/analysis/MOFA_results_plotting_and_enrichment.md)
