# Summary
## 1. Preparation
As a preparation for the data analysis workflow, thoroughly read the given manuscripts. 
- What are the advantages of the scRNA-Seq method compared to bulk RNA-Seq?
bulk RNA-Seq measures "only" the average gene expression across the population of cells in a sample. It is possible to identify differences between sample conditions. 
scRNA-Seq measures the gene expression of individual cells in a sample. It is possible to identify differences between all cell types/states. 

- What are the basic steps of the scRNA-Seq analysis?

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

- What are common problems and how are they typically solved?


- What are the major challenges in integrating single-cell transcriptomic data across different conditions, technologies, and species? How they can be solved?
- What problem does the workflow at hand address (the Seurat vignette linked above)?
## 2. Replication 
to replicate the tutorial, you need to reproduce all figures presented in the workflow. Address at least the following questions:
- Is a replication of the tutorial possible? Compare the tutorial against the rules/recommendations from Sandve et al. 2013.; comment on the clarity of the description and documentation.
- How did you set up the required environment? 
- Explain all the steps of the vignette in your own words. 
## Expanding the work
Find a publicly available data set and apply the same workflow. You may need to adapt some of the code to make it work. 
- What challenges did you ace when applying the workflow to a new data set?
- What code modifications were required?
- Are the results comparable to the results of the original tutorial, or do they deviate in some unexpected ways?
- Discuss all the results and interpret them. 

