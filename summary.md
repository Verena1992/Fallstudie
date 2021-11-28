# Summary

## 1. Preparation
As a preparation for the data analysis workflow, thoroughly read the given manuscripts. 

### What are the advantages of the scRNA-Seq method compared to bulk RNA-Seq?
Bulk RNA-Seq measures "only" the average gene expression across the population of cells in a sample. It is possible to identify differences between sample conditions. 
ScRNA-Seq measures the gene expression of individual cells in a sample.
Sequencing of individual cells simply provides the opportunity to view gene expression at a very high resolution.
It is possible to identify differences between all cell types/states. 

### What are the basic steps of the scRNA-Seq analysis?
###Should we shortly mention the steps of the experimental workflow?
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
In general, the unification and standardisation of a workflow. The choice of the programming language is also often a choice between analysis tools.
There are now very popular platforms such as Seurat, Scater or Scanpy which offer integrated working environments with large toolboxes.
In experimental workflow there are, of course, commonly known "problems" such as multiple cells beeing captured together, so doublets or multiplets can arise or that no cells are captured at all.
We mention this because it can lead to misinterpretation in the course of the data analysis. For example, unexpected high counts and a huge number of detected genes may represent doublets.
It is always an advantage to be well informed about the experimental approach, which often facilitates data interpretation.

Quality control:
 - 3 QC covariates:
   - Number of counts per barcode (count depth)
   - Number of genes per barcode
   - Fraction of counts from mitochondrial genes per barcode
Considering one of this QC covariates alone/on its own could lead to misinterpretation of cellular signals and setting the "wrong" treshold.
Therefore, these parameters should always be considered together and thresholds should be set as permissive as possible to avoid filtering out viable cell populations.
Future approach - Filtering models that account for multivariate QC dependencies provide more sensitive options in future.
What should always be taken into account is the heterogeneity in biological samples.
So heterogeneous mixtures of cell types could misinterpreted as low-quality outlier cells.
Summarized, it is relevant to set the "right" threshold and it could be helpful to the effects of quality control decissions several times.
One guidline is to use the minimum cell cluster size that is of interest.
Also the threshold should scale with the number of cells in the dataset and the intended downstream analysis.

When comparing gene expression between cells based on count data, a difference may have arisen due to sampling effects.
Therefore normalization is performed. By scaling the count data and thus determining correkt relative gene expression abundances between cells.
The most used is count depth scaling, also referred to as "counts per million" or "CPM normalization".
Normalized data should be log(x+1)-transformed for use with downstream analysis methods that assume data are normally distributed.

Batch correction:
A recommendation performing batch correction is using ComBat when cell type and state compositions between batches are consistent.
ComBat consists of a linear model of gene expression. Here, the contribution of the batches is taken into account in the mean and the variance of the data.
Data integration and batch correction should be performed by different methods.

Visualization and summarization:
Dimensionality reduction methods should be considered separately for visualization and summarization.
For exploratory visualization UMAP is a recommondation.
PCA in general for summarization  or diffusion maps for trajectory inference summarization.


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
The source of the tutorial is a script stored in a GitHub repository, so they used Version Control (Rule 4). The script reads in the raw data directly, data is not modified manualy (Rule 2). 

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
- Setup the Seurat Object
  The transcriptome of 2700 cells was measured with Illumina NextSeq 500,
  matrix.mtx from the raw Data:
  
  %%MatrixMarket matrix coordinate real general
  %
32738 2700 2286884
32709 1 4
32707 1 1
32706 1 10
32704 1 1
  
<<<<<<< HEAD
  First line shows nfeatures = 32728, ncells=2700 and number of features measured all the cells = 2286884
  Second line represent that from feature32709 4 molecules was detected in cell1  
=======
  first line shows nfeatures = 32728, ncells=2700 and number of features measured in all the cells = 2286884
  second line represent that from feature32709 4 molecules was detected in cell1  
>>>>>>> 5fa31eb759a54c7a726d6dbea7eeffe9a8c4b7c7
  
  With the Read10X() function it is possible to transform the raw data into a matrix (ncells = ncolms, ngenes = nrow), here we get 2700 columns and 32738 rows (88392600 elements). With  CreateSeuratObject() function this matrix in than stored in a Seurat object, in a Seurat object it is possible to add to the matrix also results of analysis. Only a subset of the original matrix is stored, features detected in less than 3 cells and cells with less than 200 features are excluded. The dimension of the new matrix = 13714 x 2700. 
  to save memory the numerous 0 are replaced by .
  
  
- Standard pre-processing workflow
  - QC and selecting cells for further analysis
  A barcode may mistakenly tag multiple cells (doublet) or may not tag any cells (empty droplet/well). 
  Barcodes with very few genes/number of molecules may not correspond to cells with good quality -- or quiescent
  Barcodes with high counts my correspond to douplets -- or larger in size 
  Barcodes with high fraction of mitochondrial counts are cells whose cytoplasmic mRNA has leaked out throught a broken membran -- or involved in respiratory processes
  Before analysing the single-cell gene expression data, cell with genes over 2,500 or less than 200 are excluded. Additionally if more than 5 % of this genes are mitochondrial genes than they are also excluded. Dimension of the matrix for further analysis = 13714 x 2638
  
  

- Normalizing the data
Count depths for identica cells can differ, therefore to compare cells normalization is needed. 
Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed

- Identification of highly variable features (feature selection)
Many genes will not be informative for clustering, to reduce dimensionality of the datasets only highly variable genes are selected. From the 13714 genes the 2000 most varable genes are stored under @var.features
- Scaling the data
To improve comparison between genes, gene counts are scaled. This scalinag has the effect that all genes are weighted equally (with this step biological information could be lost)

- Perform linear dimensional reduction
Biological information in expression data can be described by far fewer dimensions than the number of genes. For that a PCA is performed with the 2000 variable features 
- Determine the ‘dimensionality’ of the dataset
With PCA it is possible to summarize a dataset. How many pc are neede to have enough information can be derminated by "elbow" heuristics, or the permutation-test-based jackstraw method.
- Cluster the cells
  Euclidean distances form cells on the PC-reduced expression space are calculated. 
  *K-Nearest Neighbour graph*
  Based on this distance, each cell (=node) is connected to 20 other cells (with the lowest distanceses). This results for 2638 cells to 52760 connections (pbmc@graphs[["RNA_nn"]]@x). 
  *Shared Nearest-neighbor graph (snn)*
  Additionaly, the neighborhood overlap (Jaccard index) between every cell and its k.param nearest neighbors is calculated. 
  
  
  
- Run non-linear dimensional reduction (UMAP/tSNE)
- Finding differentially expressed features (cluster biomarkers)
- Assigning cell type identity to clusters

## Expanding the work
Find a publicly available data set and apply the same workflow. You may need to adapt some of the code to make it work. 
- What challenges did you ace when applying the workflow to a new data set?
- What code modifications were required?
- Are the results comparable to the results of the original tutorial, or do they deviate in some unexpected ways?
- Discuss all the results and interpret them. 

