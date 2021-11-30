# Summary

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
   - 
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
  
  Possible changes in cellular density (shifts in subpopulation frequency) between conditions, 
  change in feature scale across conditions (global transcriptional shifts, differences in normalization strategies)
  non-overlapping populations
  
  For identification of shared population across data sets, subpopulations must be aligned. With canonical correlation analysis, shared correlation structures across data sets are identifyed. Basic vectors from CCA are linear transformed to correct global shifts or normalization strategy. To correct changes in population density (nonlinear shift) dynamic time warping is used.
  
  To identify non-overlapping populations: Because CCA may struggle to identify rare subpopulations present in only one data set, PCA  may be able to separate these cells. PCA is performed on each data set independently and explained variance compared with CCA.

### What problem does the workflow at hand address (the Seurat vignette linked above)?
Clustering - Marker identification - Cluster annotation

## 2. Replication 
To replicate the tutorial, you need to reproduce all figures presented in the workflow. Address at least the following questions:

### Is a replication of the tutorial possible? Compare the tutorial against the rules/recommendations from Sandve et al. 2013.; comment on the clarity of the description and documentation.



In general, the tutorial and the results seems reproducible, but in the jackstraw plot we get different numerical values for the PC7, PC9 and PC12. Here, however, the question is whether this is relevant at all, since no significant changes are involved.
The course of the line do not change and they are all above the dashed line.
Original tutorial: PC7 3.38e-23, PC9 5.22e-12, PC12 5.92e-05
Reproduced tutorial: PC7 4.05e-24, PC9 8.8e-12, PC12 0.000101
When displaying/arranging the heat maps of the 15 main components. We also noticed that the arrangement of the features/genes in the tutorial does not correspond to the reproduced one. See PC 15, XCL2 last feature according to the tutorial. In the reproduced tutorial APOBEC3H is the last feature.
What we also noticed is that the output of the pvalues is different, although we actually use the exact code from the tutorial. 
In the seurat repository we found a discussion with exactly this issue [seurat/issue1789](https://github.com/satijalab/seurat/issues/1789), mentioned that "It is possible for these values to change slightly, based on the exact computing platform." And that means it is not 100 % reproducibel. 

As well in the original script there is a little bug, for the chunk calculating the object size, the chunk is made with ''''r instad of '''r resulting in a error message when trying to convert into html. 


To ensure our results, we ran our reproduction on 2 different computers and executed the commands in the tutorial individually as well as in a complete run. No changes were apparent to us here.

- Rule 1: For every result, keep track of how it was produced
  - We would have done it that way and it was obviously feasible for us
- Rule 2: Avoid manual data manipulation steps
  - No manual data manipulation was carried out, script reads in the raw data directly
- Rule 3: Archive the exact version of all external programs used
 with Session Info the exact names and versions of the main programs are given, and they also provide a docker images. Instructions on how to set up required environments or how to use docker is missing. 
used IDE is also mentioned. 
- Rule 4: Version Control All Custom Scripts
  - The source of the tutorial is a script stored in a GitHub repository, so they used Version Control
- Rule 5: Record all intermediate results, when possible in standardized formats
  - at the end "saveRDS(pbmc, file = "../data/pbmc3k_final.rds")" the seurat object is stored. But the object serves as a container with different slot for data and analysis, that makes it possible to see for example scaled and not scaled data without overwriting. 
- Rule 6: For analyses that include randomness, note underlying random seeds
  - for quiet a lot of functions, random.seeds are set, but newer mentioned directly. For the function JackStraw(pbmc, num.replicate = 100) randomness is mentioned but, repeated execution doesn not change the result. In the source code we found that also here a random.seed is set. 
- Rule 7: Always store raw data behind plots
  - The corresponding data is stored
- Rule 8: Generate hierarchical analysis output, allowing layers of increasing detail to be inspected
  - we found that the seuratObject, is a format, that allows to inspect different layers of the result. 
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
  
first line shows nfeatures = 32728, ncells=2700 and number of features measured in all the cells = 2286884 \
second line represent that from feature32709 4 molecules was detected in cell1 

With the Read10X() function it is possible to transform the raw data into a matrix (ncells = ncolms, ngenes = nrow), here we get 2700 columns and 32738 rows (88392600 elements). With  CreateSeuratObject() function this matrix in than stored in a Seurat object, in a Seurat object it is possible to add to the matrix also results of analysis. Only a subset of the original matrix is stored, features detected in less than 3 cells and cells with less than 200 features are excluded. The dimension of the new matrix = 13714 x 2700. 
To save memory the numerous 0 are replaced by .(point)
  
- **QC and selecting cells for further analysis**

 A barcode may mistakenly tag multiple cells (doublet) or may not tag any cells (empty droplet/well). \
 Barcodes with very few genes/number of molecules may not correspond to cells with good quality -- or quiescent \
 Barcodes with high counts my correspond to douplets -- or larger in size \
 Barcodes with high fraction of mitochondrial counts are cells whose cytoplasmic mRNA has leaked out throught a broken membran -- or involved in respiratory processes \
Before analysing the single-cell gene expression data, cell with genes over 2,500 or less than 200 are excluded. Additionally if more than 5 % of this genes are mitochondrial genes than they are also excluded. Dimension of the matrix for further analysis = 13714 x 2638
  
- **Normalizing the data**

Count depths for identical cells can differ, therefore to compare cells normalization is needed. 
Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed

- **Identification of highly variable features (feature selection)**

Many genes will not be informative for clustering, to reduce dimensionality of the datasets only highly variable genes are selected. From the 13714 genes the 2000 most varable genes are stored under @var.features

- **Scaling the data**

To improve comparison between genes, gene counts are scaled. This scalinag has the effect that all genes are weighted equally (with this step biological information could be lost)

- **Perform linear dimensional reduction**

Biological information in expression data can be described by far fewer dimensions than the number of genes. For that a PCA is performed with the 2000 variable features 

- **Determine the ‘dimensionality’ of the dataset**

With PCA it is possible to summarize a dataset. How many pc are neede to have enough information can be derminated by "elbow" heuristics, or the permutation-test-based jackstraw method.

- **Cluster the cells**

Euclidean distances form cells on the PC-reduced expression space are calculated. 

*K-Nearest Neighbour graph* \
Based on this distance, each cell (=node) is connected to 20 other cells (with the lowest distanceses). This results for 2638 cells to 52760 connections (pbmc@graphs[["RNA_nn"]]@x). 

*Shared Nearest-neighbor graph (snn)* \
Additionaly, the neighborhood overlap (Jaccard index) between every cell and its k.param nearest neighbors is calculated. 
  
Then with the SNN graph, the modularity function is optimized to determine the clusters. For that the Louvain algorithm is used. Graphs with a high modularity score will have many connections within a community but only few pointing outwards to other communities. The Louvain algorithm detects communities as groups of cells that have more links between them than expected from the number of links the cells habe in total, or in other words: the fraction of edges that run within the community is compared to the fraction we would expect to find if we "randomly rewired" the network. 
  
  
  
- **Run non-linear dimensional reduction (UMAP/tSNE)**

cells within the same cluster should group together 

in the tutorial not the python package umap is used, instead the r package uwot is used. Installation of python package is not needed.

> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric /
> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'

- **Finding differentially expressed features (cluster biomarkers)**

To detected genes that are expressed much more/less in one cluster compared to the other cells.

ident.1 = cluster to compare (if ident.2 is not defined, all other cells are used for comparison)
idnet.2 = A second identity class for comparison. ident.1 vs ident.2 

with the function FindAllMarkers(), markers for all cluster are determined. 

if limma is not installed in environment, it is installed automatically with:

> install.packages('BiocManager')
> BiocManager::install('limma')

to speed up this process there are 3 options:
1. set the min.pct argument to a higher value than the default of 0.1 (excludes genes that are very infrequnetly expressed)
2. set the logfc.threshold argument to a higher value than the default of 0.25 (excludes weaker signals)
3. set the max.cells.per. ident argument to a number, default is Inf. (Down sample each cluster to a max number of cells)


to calculate the impact a marker gene has on the classification the ROC test can be used.



- **Assigning cell type identity to clusters**

markers are used to define the cell type of the clusters. 

## Expanding the work
Find a publicly available data set and apply the same workflow. You may need to adapt some of the code to make it work.
We used the Human Glioblastoma Multiforme:3'v3 Targeted, Neuroscience Panel
Single Cell Gene Expression Dataset by Cell Ranger 4.0.0
https://www.10xgenomics.com/resources/datasets/human-glioblastoma-multiforme-3-v-3-targeted-neuroscience-panel-3-standard-4-0-0
Human Glioblastoma Multiforme cells from a male donor aged 57 were obtained by 10x Genomics from Discovery Life Sciences.
The Sequencing was done by Illumina NovaSeq 6000
filtered 1138 features across 4433 samples within 1 assay
raw 5277 features across 4616 samples within 1 assay

- **What challenges did you face when applying the workflow to a new data set?**
  - First you have to think about which biological data do you use and check the raw data preparation (like single cell sequencing or bulk sequencing)

  - Make sure that you use the correct dataset (raw data or filtered dataset)
    The cellranger pipeline outputs two types of feature-barcode matrices:
    Unfiltered feature-barcode matrix and Filtered feature-barcode matrix.
    The filtered matrix contains only detected cellular barcodes. For Targeted Gene Expression samples, non-targeted genes are removed from the filtered matrix).
	Because it was shown that using targeted gene expression with the Human Neuroscience Panel (Number of Genes Targeted = 1186) enables comprehensive and efficient characterization of the human brain and nervous system, we decided to work with targeted data. 
	All analysis for targeted data is performed with the filtered matrix (with removed non-targeted genes).
  - Check which input files and fileformats are necessary for the current workflow/pipeline
  - Make sure that you do not mix or transfer the already used dataset with the new one
- finding a suited parameter (min.cells, min.features) for creating Seurat Object as measured features of the cell is lower than in the tutorial
  - with the targeted approach no genes starting with 'MT-' are present in the dataset. So mitochondiral contamination as signal of low-quality cells is not available. 
- because the target approach already tries to focus on genes that matter most, feature number is low compared to the tutorial. The step feature selection is may not needed. 
- it seems that the dimensionality is higher
- Assigning cell type identity to clusters not possible without having canonical markers





- **What code modifications were required?**
  - name from the read in data
  - lower number of min.features for cells beeing included in seurat Object 
  - You have to adapt some filters, because you analyze different number of features and samples
  - increase number of used dimenstion
  - skip feature selection step
  - because high heterogenity is expected also resolution parameter is set to a higher number
  - change the variable name pbmc, because data is from different cell type. To make it more robust for different datasets a neutral name (scSEQ) is used. 

- **Are the results comparable to the results of the original tutorial, or do they deviate in some unexpected ways?**
  - The results cannot be compared in the biological context because the raw data are completely different and, above all, the cells are completely different. Tumour cells are basically very heterogeneous and in the case of glioblastoma cells this applies even more extensively.
 - because a panel was used, in the matrix also features (n=20) were present without beeing detected in a cell. 
  - The many black areas in the heat maps are not yet clear to me???? low = "magenta", high = "yellow", mid = "black" 
  - In our data set, there are also apparently more principal components beforehand due to the heterogeneity of the cells, and the individual PCs are not so easily separable.
  - Whereby the "main" elbow again at 8 - 9 PCs
  - But if you compare the plots to the basic tutorial, there are similarities.
  - The QC violin plots with the mitochondrial RNA look the most different.

- **Discuss all the results and interpret them.**
  - Basically, it is not so easy to interpret the available data because we have not done it to this extent before. Above all, there is the question of clinical relevance.
  - At first glance, the up- and down-regulated genes in the specific cell variants can be viewed.


