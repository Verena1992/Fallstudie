# Summary
## 1. Preparation
As a preparation for the data analysis workflow, thoroughly read the given manuscripts. 
### What are the advantages of the scRNA-Seq method compared to bulk RNA-Seq?
bulk RNA-Seq measures "only" the average gene expression across the population of cells in a sample. It is possible to identify differences between sample conditions. 
scRNA-Seq measures the gene expression of individual cells in a sample. It is possible to identify differences between all cell types/states. 

### What are the basic steps of the scRNA-Seq analysis?

- pre-processing:
  - quality control
  - normalization
  - data correction and integration
  - feature selection, dimensionality reduction and visualization
- downstream analysis:
  - cell level and gene-level
    - cluster analysis
     - clustering (cell-level)
     - cluster annotation (gene-level)
     - compositional analysis (cell-level)
    - trajectory analysis
     - trajectory inference (cell-level)
     - gene expression dynamics (gene-level)
     - metastable states (cell-level)
 
  - gene level
    - differential expression analysis
    - gene set analysis
    - gene regulatory networks

### What are common problems and how are they typically solved?


### What are the major challenges in integrating single-cell transcriptomic data across different conditions, technologies, and species? How they can be solved?
  it is difficult to distinguish between changes in the proportional composition of cell types in a sample and expression changes within a given cell type.
  
  possible changes in cellular density (shifts in subpopulation frequency) between conditions, 
  change in feature scale across conditions (global transcriptional shifts, differences in normalization strategies)
  non-overlapping populations
  
  for identification of shared population across data sets, subpopulations must be aligned. With canonical correlation analysis, shared correlation structures across data sets are identifyed. Basic vectors from CCA are linear transformed to correct global shifts or normalization strategy. To correct changes in population density (nonlinear shift) dynamic time warping is used.
  
 to identify non-overlapping populations: because cca may struggle to identify rare subpopulations present in only one data set, pca may be able to separate these cells.Pca is performed on each data set independently and explained variance compared with cca.
### What problem does the workflow at hand address (the Seurat vignette linked above)?
Clustering - Marker identification - Cluster annotation
## 2. Replication 
to replicate the tutorial, you need to reproduce all figures presented in the workflow. Address at least the following questions:
### Is a replication of the tutorial possible? Compare the tutorial against the rules/recommendations from Sandve et al. 2013.; comment on the clarity of the description and documentation.
the source of the tutorial is a script stored in a GitHub repository, so they used Version Control (Rule 4). The script reads in the raw data directly, data is not modified manualy (Rule 2). 
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
  the transcriptome of 2700 cells was measured with illumina,
  matrix.mtx from the raw Data:
  
  %%MatrixMarket matrix coordinate real general
  %
32738 2700 2286884
32709 1 4
32707 1 1
32706 1 10
32704 1 1
  
  first line shows nfeatures = 32728, ncells=2700 and number of features measured all the cells = 2286884
  second line represent that from feature32709 4 molecules was detected in cell1  
  
  with the Read10X() function it is possible to transform the raw data into a matrix (ncells = ncolms, ngenes = nrow), here we get 2700 columns and 32738 rows (88392600 elements). With  CreateSeuratObject() function this matrix in than stored in a Seurat object, in a Seurat object it is possible to add to the matrix also results of analysis. Only a subset of the original matrix is stored, features detected in less than 3 cells and cells with less than 200 features are excluded. The dimension of the new matrix = 13714 x 2700. 
  to save memory the numerous 0 are replaced by .
  
  
- Standard pre-processing workflow
  - QC and selecting cells for further analysis
  A barcode may mistakenly tag multiple cells (doublet) or may not tag any cells (empty droplet/well). 
  Barcodes with very few genes/number of molecules may not correspond to cells with good quality -- or quiescent
  Barcodes with high counts my correspond to douplets -- or larger in size 
  Barcodes with high fraction of mitochondrial counts are cells whose cytoplasmic mRNA has leaked out throught a broken membran -- or involved in respiratory processes
  Before analysing the single-cell gene expression data, cell with genes over 2,500 or less than 200 are excluded. Additionally if more than 5 % of this genes are mitochondrial genes than they are also excluded. Dimension of the matrix for further analysis = 13714 x 2638
  
  

- Normalizing the data
count depths for identica cells can differ, therefore to compare cells normalization is needed. 
Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed

- Identification of highly variable features (feature selection)
many genes will not be informative for clustering, to reduce dimensionality of the datasets only highly variable genes are selected. From the 13714 genes the 2000 most varable genes are stored under @var.features
- Scaling the data
To improve comparison between genes, gene counts are scaled. This scalinag has the effect that all genes are weighted equally (with this step biological information could be lost)

- Perform linear dimensional reduction
biological information in expression data can be described by far fewer dimensions than the number of genes. For that a PCA is performed with the 2000 variable features 
- Determine the ‘dimensionality’ of the dataset
with PCA it is possible to "summarize" a dataset. How many pc are neede to have enough information can be derminated by "elbow" heuristics, or the permutation-test-based jackstraw method.
- Cluster the cells
- Run non-linear dimensional reduction (UMAP/tSNE)
- Finding differentially expressed features (cluster biomarkers)
- Assigning cell type identity to clusters

## Expanding the work
Find a publicly available data set and apply the same workflow. You may need to adapt some of the code to make it work. 
- What challenges did you ace when applying the workflow to a new data set?
- What code modifications were required?
- Are the results comparable to the results of the original tutorial, or do they deviate in some unexpected ways?
- Discuss all the results and interpret them. 

