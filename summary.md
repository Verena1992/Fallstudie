# Summary
- Gebhart Verena, 118787
- Hennebichler Bernhard, 116252
- You can find our github repository here [Fallstudie](https://github.com/Verena1992/Fallstudie)

## 1. Preparation
As a preparation for the data analysis workflow, thoroughly read the given manuscripts. 

### What are the advantages of the scRNA-Seq method compared to bulk RNA-Seq?
Bulk RNA-Seq measures "only" the average gene expression across the population of cells in a sample. It is possible to identify differences between sample conditions. 
ScRNA-Seq measures the gene expression of individual cells in a sample.
Sequencing of individual cells simply provides the opportunity to view gene expression at a very high resolution.
It is possible to identify differences between all cell types/states. 

### What are the basic steps of the scRNA-Seq analysis?
- Experimental workflow before analysis:
  - Single-cell dissociation
  - Single-cell isolation
  - Library construction
  - Sequencing
 
- Pre-processing:
  - Quality control
  - Normalization
  - Data correction and integration
  - Feature selection, dimensionality reduction and visualization

- Downstream analysis:
 - Cell-level and gene-level
    - Cluster analysis
      - Clustering (cell-level)
      - Cluster annotation (gene-level)
      - Compositional analysis (cell-level)
    - Trajectory analysis
      - Trajectory inference (cell-level)
      - Gene expression dynamics (gene-level)
      - Metastable states (cell-level)
 
 - Gene-level
    - Differential expression analysis
    - Gene set analysis
    - Gene regulatory networks

### What are common problems and how are they typically solved?
First, the unification and standardisation of an analysis workflow. The choice of the programming language is often a choice between analysis tools.
There are now very popular platforms such as Seurat, Scater or Scanpy which offer integrated working environments with large toolboxes.
In experimental workflow there are commonly known "problems" such as multiple cells beeing captured together. So doublets or multiplets can arise or no cells are captured at all. This could lead to misinterpretation in the course of the data analysis. For example, unexpected high counts and a huge amount of detected genes may represent doublets.
It is benefical to be well informed about the experimental approach/raw data, which often facilitates data interpretation.

Quality control - 3 QC covariates:
   - Number of counts per barcode (count depth)
   - Number of genes per barcode
   - Fraction of counts from mitochondrial genes per barcode
   
Considering one of this QC covariates alone could lead to misinterpretation of cellular signals and setting the "wrong" threshold.
Therefore, these parameters should always be considered together and thresholds should be set as permissive as possible to avoid filtering out viable cell populations.
Future approach - Filtering models that account for multivariate QC dependencies provide more sensitive options.
The heterogeneity of biological samples should always be taken into account in data analysis.
One example - Heterogeneous mixtures of cell types could misinterpreted as low-quality outlier cells.
Summarized, it is relevant to set the threshold carefully. It could be helpful to check the effects of quality control decissions several times.
One guideline is to use the minimum cell cluster size that is of interest.
Also the threshold should scale with the number of cells in the dataset and the intended downstream analysis.

When comparing gene expression between cells based on count data, a difference may have arisen due to sampling effects.
Therefore normalization is performed. By scaling the count data and thus determining correct relative gene expression abundances between cells.

Batch correction:
A recommendation performing batch correction is using ComBat when cell type and state compositions between batches are consistent.
ComBat consists of a linear model of gene expression. Here, the contribution of the batches is taken into account in the mean and the variance of the data.
Data integration and batch correction should be performed by different methods.

Visualization and summarization:
Dimensionality reduction methods should be considered separately for visualization and summarization.
For exploratory visualization UMAP is a recommondation.
PCA in general for summarization or diffusion maps for trajectory inference summarization.


### What are the major challenges in integrating single-cell transcriptomic data across different conditions, technologies, and species? How they can be solved?
It is difficult to distinguish between changes in the proportional composition of cell types in a sample and expression changes within a given cell type.
  
 - Possible changes in cellular density (shifts in subpopulation frequency) between conditions 
 - Change in feature scale across conditions (global transcriptional shifts, differences in normalization strategies)
 - Non-overlapping populations
  
For identification of shared population across data sets, subpopulations must be aligned. With canonical correlation analysis, shared correlation structures across data sets are identifyed. Basic vectors from CCA are linear transformed to correct global shifts or normalization strategy. To correct changes in population density (nonlinear shift) dynamic time warping is used.
  
To identify non-overlapping populations: Because CCA may struggle to identify rare subpopulations present in only one data set, PCA  may be able to separate these cells. PCA is performed on each data set independently and explained variance compared with CCA.

### What problem does the workflow at hand address (the Seurat vignette linked above)?
Clustering --> Marker identification --> Cluster annotation

## 2. Replication 
To replicate the tutorial, you need to reproduce all figures presented in the workflow. Address at least the following questions:

### Is a replication of the tutorial possible? Compare the tutorial against the rules/recommendations from Sandve et al. 2013.; comment on the clarity of the description and documentation.

In general, the tutorial and the results seems reproducible, but in the jackstraw plot we get different numerical values for PC7, PC9 and PC12.
Original tutorial: PC7 3.38e-23, PC9 5.22e-12, PC12 5.92e-05
Reproduced tutorial: PC7 4.05e-24, PC9 8.8e-12, PC12 0.000101
However no significant changes are involved. The course of the line do not change obviously and they are all above the dashed line.
When displaying/arranging the heat maps of the 15 main components. We also noticed that the arrangement of the features/genes in the tutorial does not correspond to the reproduced one. For example PC 15, XCL2 last feature according to the tutorial. In the reproduced tutorial APOBEC3H is the last feature.
What we also noticed is that the output of the pvalues is different, although we actually use the exact code from the tutorial. 
In the seurat repository we found a discussion with exactly this issue [seurat/issue1789](https://github.com/satijalab/seurat/issues/1789), mentioned that "It is possible for these values to change slightly, based on the exact computing platform." And that means it is not 100 % reproducible. 

As well in the original script there is a little bug, for the chunk calculating the object size, the chunk is made with ''''r instad of '''r resulting in a error message when trying to convert into html. 

To ensure our results, we ran our reproduction on 2 different computers and executed the commands in the tutorial individually as well as in a complete run. No changes were apparent to us here.

- Rule 1: For every result, keep track of how it was produced
  - We would have done it that way and it was obviously feasible for us
- Rule 2: Avoid manual data manipulation steps
  - No manual data manipulation was carried out, script reads in the raw data directly
- Rule 3: Archive the exact version of all external programs used
  - With Session Info the exact names and versions of the main programs are given, and they also provide a docker image. Instructions on how to set up required environments or how to use docker is missing. Used IDE is also mentioned. 
- Rule 4: Version Control All Custom Scripts
  - The source of the tutorial is a script stored in a GitHub repository, so they used Version Control
- Rule 5: Record all intermediate results, when possible in standardized formats
  - At the end "saveRDS(pbmc, file = "../data/pbmc3k_final.rds")" the seurat object is stored. But the object serves as a container with different slot for data and analysis, that makes it possible to see for example scaled and not scaled data without overwritting. 
- Rule 6: For analyses that include randomness, note underlying random seeds
  - For quiet a lot of functions, random.seeds are set, but never mentioned directly. For the function JackStraw(pbmc, num.replicate = 100) randomness is mentioned but, repeated execution does not change the result. In the source code we found that also here a random.seed is set. 
- Rule 7: Always store raw data behind plots
  - The corresponding data is stored
- Rule 8: Generate hierarchical analysis output, allowing layers of increasing detail to be inspected
  - We found that the seuratObject, is a format, that allows to inspect different layers of the result. 
- Rule 9: Connect textual statements to underlying results
  - The tutorial is commented out in detail and linked to the necessary sources
- Rule 10: Provide public access to scripts, runs and results
  - Via Github

### How did you set up the required environment? 
**conda**

```bash
conda --version # conda 4.10.3
```

an environment with seurat v4 was created:

```bash
conda create --name seurat
conda activate seurat
conda install r-base=4.1.0 
conda install r-ggplot2=3.3.5 r-patchwork=1.1.1 r-seuratObject=4.0.2 r-seurat=4.0.4 r-dplyr=1.0.7
conda install r-rmarkdown r-knitr r-xfun r-htmltools
conda install -c bioconda bioconductor-limma=3.48.0
```

and the environment exported to a YAML file:

```bash
conda env export --name seurat > seurat.yml 
```
**docker**

*Installation*

```bash
sudo apt-get update
```

```bash
sudo apt-get remove docker docker-engine docker.io
```

```bash
sudo apt install docker.io
```

```bash
sudo systemctl start docker
sudo systemctl enable docker
```

```bash
docker --version #Docker version 20.10.7, build 20.10.7-0ubuntu5~20.04.2
```


*Build an image from a Dockerfile*

```bash
cd <path to Dockerfile>
```

```bash
sudo docker build -f Dockerfile -t seurat .
```


### Explain all the steps of the vignette in your own words. 
- **Setup the Seurat Object** \
The transcriptome of 2700 cells was measured with Illumina NextSeq 500, 
  
 matrix.mtx from the raw Data:
  
> %%MatrixMarket matrix coordinate real general \
> % \
> 32738 2700 2286884 \
> 32709 1 4 \
> 32707 1 1 \
> 32706 1 10 \
> 32704 1 1 
  
First line shows nfeatures = 32728, ncells=2700 and number of features measured in all the cells = 2286884 \
Second line represent that from feature32709 4 molecules was detected in cell1 

With the Read10X() function it is possible to transform the raw data into a matrix (ncells = ncolms, ngenes = nrow), here we get 2700 columns and 32738 rows (88392600 elements). With  CreateSeuratObject() function the matrix in stored in a Seurat object, here it is possible to add results of the analysis to the matrix. Only a subset of the original matrix is stored, features detected in less than 3 cells and cells with less than 200 features are excluded. The dimension of the new matrix = 13714 x 2700. 
To save memory the numerous 0 are replaced by .(point)
  
- **QC and selecting cells for further analysis**

 A barcode may mistakenly tag multiple cells (doublet) or may not tag any cells (empty droplet/well). \
 Barcodes with very few genes/number of molecules may not correspond to cells with good quality -- or quiescent \
 Barcodes with high counts may correspond to douplets -- or larger in size \
 Barcodes with high fraction of mitochondrial counts are cells whose cytoplasmic mRNA has leaked out throught a broken membran -- or involved in respiratory processes \
Before analysing the single-cell gene expression data, cell with genes over 2,500 or less than 200 are excluded. Additionally if more than 5 % of this genes are mitochondrial genes than they are also excluded. Dimension of the matrix for further analysis = 13714 x 2638
  
- **Normalizing the data**

Count depths for identical cells can differ, therefore to compare cells normalization is needed. 
Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed

- **Identification of highly variable features (feature selection)**

Many genes will not be informative for clustering. To reduce dimensionality of the datasets only highly variable genes are selected. From the 13714 genes the 2000 most variable genes are stored under @var.features

- **Scaling the data**

To improve comparison between genes, gene counts are scaled. This scaling has the effect that all genes are weighted equally. With this step biological information could be lost.

- **Perform linear dimensional reduction**

Biological information in expression data can be described by far fewer dimensions than the number of genes. For that a PCA is performed with the 2000 variable features. 

- **Determine the ???dimensionality??? of the dataset**

With PCA it is possible to summarize a dataset, as well as to check how many PCs are needed to obtain enough information. This can be determined by "elbow" heuristics, or the permutation-test-based jackstraw method.

- **Cluster the cells**

Euclidean distances form cells on the PC-reduced expression space are calculated. 

*K-Nearest Neighbour graph* \
Based on this distance, each cell (=node) is connected to 20 other cells (with the lowest distanceses). This results for 2638 cells to 52760 connections (pbmc@graphs[["RNA_nn"]]@x). 

*Shared Nearest-neighbor graph (snn)* \
Additionaly, the neighborhood overlap (Jaccard index) between every cell and its k.param nearest neighbours is calculated. 
  
Then with the SNN graph, the modularity function is optimized to determine the clusters. For that the Louvain algorithm is used. Graphs with a high modularity score will have many connections within a community but only few pointing outwards to other communities. The Louvain algorithm detects communities as groups of cells that have more links between them than expected from the number of links the cells have in total, or in other words: The fraction of edges that run within the community is compared to the fraction we would expect to find if we "randomly rewired" the network. 
  
  
- **Run non-linear dimensional reduction (UMAP/tSNE)**

Cells within the same cluster should group together. 

The tutorial do not use the python package umap, instead the r package uwot is used. Installation of python package is not needed.

> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric /
> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'

- **Finding differentially expressed features (cluster biomarkers)**

To detect genes that are expressed much more or less in one cluster compared to the other cells.

ident.1 = cluster to compare (if ident.2 is not defined, all other cells are used for comparison
idnet.2 = A second identity class for comparison. ident.1 vs ident.2 

With the function FindAllMarkers(), markers for all cluster are determined. 

If limma is not installed in environment, it is installed automatically with:

> install.packages('BiocManager')
> BiocManager::install('limma')

To speed up this process there are 3 options:
1. Set the min.pct argument to a higher value than the default of 0.1 (excludes genes that are very infrequently expressed)
2. Set the logfc.threshold argument to a higher value than the default of 0.25 (excludes weaker signals)
3. Set the max.cells.per. ident argument to a number, default is Inf. (Down sample each cluster to a max number of cells)

To calculate the impact of a marker gene to the classification the ROC test can be used.

- **Assigning cell type identity to clusters**

Markers are used to define the cell type of the clusters. 

## Expanding the work
Find a publicly available data set and apply the same workflow. You may need to adapt some of the code to make it work.

We used the Human Glioblastoma Multiforme:3'v3 Targeted, Neuroscience Panel
Single Cell Gene Expression Dataset by Cell Ranger 4.0.0
Provided by 10xGenomics:
https://www.10xgenomics.com/resources/datasets/human-glioblastoma-multiforme-3-v-3-targeted-neuroscience-panel-3-standard-4-0-0
Human Glioblastoma Multiforme cells from a male donor aged 57 were obtained by 10x Genomics from Discovery Life Sciences.
The Sequencing was done by Illumina NovaSeq 6000
filtered 1138 features across 4433 samples within 1 assay
raw 5277 features across 4616 samples within 1 assay

- **What challenges did you face when applying the workflow to a new data set?**

  - At the beginning it is useful to obtain information of the examined biological samples and check the raw data preparation (like single cell sequencing or bulk sequencing)
  - Make sure that you use the correct dataset (raw data or filtered dataset)
    The cellranger pipeline outputs two types of feature-barcode matrices:
    Unfiltered feature-barcode matrix and Filtered feature-barcode matrix.
    The filtered matrix contains only detected cellular barcodes. For Targeted Gene Expression samples, non-targeted genes are removed from the filtered matrix.
    It was shown that using targeted gene expression with the Human Neuroscience Panel (Number of Genes Targeted = 1186) enables comprehensive and efficient characterization of the human brain and nervous system. We decided to work with targeted data. All analysis for targeted data is performed with filtered matrix (with removed non-targeted genes).
  - Check which input files and fileformats are necessary for the current workflow/pipeline
  - Make sure that you do not mix or transfer the already used dataset with the new one
  - Finding a suited parameter (min.cells, min.features) for creating Seurat Object as measured features of the cell is lower than in the tutorial.
  - With the targeted approach no genes starting with 'MT-' are present in the dataset. So mitochondiral contamination as signal of low-quality cells is not available. 
  - The step feature selection is may not needed, because the target approach already tries to focus on genes that matter most.
  - Feature number is low compared to the tutorial, but it seems the dimensionality is higher.
  - Assigning cell type identity to clusters not possible without having canonical markers

- **What code modifications were required?**
  - Name from the read in data
  - Lower number of min.features for cells beeing included in seurat Object 
  - You have to adapt some filters, because you analyze different number of features and samples
  - Increase number of used dimension
  - Skip feature selection step
  - Because high heterogenity (due to biological context) is expected also resolution parameter is set to a higher count
  - Change the variable name pbmc, because data is from different cell type. To make it more robust/adaptable for different datasets a neutral name (scSEQ) is used. 

- **Are the results comparable to the results of the original tutorial, or do they deviate in some unexpected ways?**
  - The results are not comparable through biological context. The existing cells are completely different. Tumour cells are basically very heterogeneous and in the case of glioblastoma cells this applies even more extensively.
  - Because a panel was used, in the matrix also features (n=20) were present without beeing detected in a cell.  
  - In our data set, there are also apparently more principal components beforehand due to the heterogeneity of the cells, and the individual PCs are not so easily separable.
  - Whereby the "main" elbow again at 8 - 9 PCs
  - But if you compare the plots to the basic tutorial, there are similarities.
  - The QC violin plots with the mitochondrial RNA look the most different.

- **Discuss all the results and interpret them.**
  - The dataset only contains cells from a glioblastoma, which means this is a neuro-targeted aprroach and not a whole transcriptome approach.
  - QC and selecting cells for further analysis:
     - As mentioned in the section above, with the targeted approach no genes starting with 'MT-' are present in the dataset. So mitochondiral contamination as signal of low-quality cells is not available, because the dataset is already filtered for healthy and alive cancer cells. It should not contain dead cells. The distribution of counts and features can be seen in the violin plots. The nFeature_RNA plot shows an even distibution. A clustered occurrence is visible from 200 to 400. The nCount_RNA plot shows a huge clustered occurence at the bottom. In the FeatureScatter feature-feature relationships are visualized. As expected, there is no mtDNA present in the data. In the right plot, nFeature_RNA is plotted against nCount_RNA. In our data, a wider scatter and steeper increase can be observed compared to the tutorial. For this reason the value of 0.88 is slightly lower compared to the original tutorial.
  - Identification of highly variable features (feature selection) - 10 most highly variable genes were identified:
     -  We want to remove any technical artifacts, therefore normalized data is used. Standardized variance is plotted against average expression. We used the default settings of 2,000 features per dataset. Therefore we have a variable count of 1166 and no non-variable count.
     -  The 10 highly variable genes are: COL4A1, APOD, MT1H, PLA2G2A, CLDN5, COL4A2, MMP9, EGFL7, CXCL10, ESAM.
  - Perform linear dimensional reduction
     - Some ways are provided to visualize cells and features that define the PCA. When PC_1 is plotted against PC_2, at first glance there are 3 clusters of glioblastoma cells with cells in between. The first two PCs explain the largest percentage of the total variance.
     - Afterwards the PCs/Dimensions were visualized with heatmaps. You can see the different gene expression levels; magenta = low gene expression , black = normal gene expression, yellow = high gene expression. If the PC clearly divides the cells into distinct groups, the borders between the clusters are clearly defined by the colors. One singal to cut off higher dimensions, when the borders start to become blurry. Compared to the tutorial results you can obsere, our dataset has also in PC13 to PC15 a variation between high and low gene expressesion, whereas in the tutorial there only low gene expressions in the last 3 visualized PCs. Which in turn confirms that in our case we are dealing with a very heterogenous/high dimensional dataset.
  - JackStrawplot and Elbow heuristics
     - JackStrawplot compares the distribution of the PCs p-value to a uniform distribution (dashed line). Significant PCs will show a strong enrichment of features with low p-values, that is strongly skewed to the left compared to the null distribution (dashed line). The p-value for each PC is based on a proportion test comparing the number of genes with a p-value below a patricular threshold, compared with the proportion of genes expected under a uniform distribution of p-values.
     - Elbowplot - Typically, the first few PCs tend to explain quite a proportion of the variance. We can observe the "main" elbow around PC 8-9, but our example isn't that clear. We can also observe smaller steps for PC 12-13 and PC 20-21.
  - Run non-linear dimensional reduction (UMAP) - Place similar cells together in low-dimensional space
     - Cells within the graph-based cluster determined above should co-localize on the dimension reduction plots. 14 cell clusters are visible.
  - Finding differentially expressed features (cluster biomarkers)
     - It identifies positive and negative markers of a single cluster
     - Visualize marker expression via ViolinPlots, we used PLP1 and MAG, the expression level is plotted against the 14 PCs.
     - Afterwards we plot the raw counts for OLR1 und PLAUR
     - Then the featureplot for the top 10 is used
     - An expression heatmap for given cells and features is created. The top 20 markers or all markers are used if less than for each cluster
  - Assigning cell type identity to cluster fortunately in the case of this dataset, we can use canonical markers to easily math the unbiased clustering known cell types.
     - We are not able to determine the canonical markers to assign cell type identity

Summarized, it can be said interpreting the data and assigning the corresponding clusters to the associated markers turned out to be more difficult than expected. In other words, our data deviates from the dataset of the original tutorial.

