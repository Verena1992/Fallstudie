# Fallstudie

an environment with seurat v4 was created:

conda create --name seurat -c conda-forge -c bioconda r-seurat=4*

and the environment exported to a YAML file:

conda activate seurat
conda env export --name seurat  > seurat.yml
conda deactivate



## Get started

### Create an environment from YAML file:

```bash
conda env create --file seurat.yml
```

```bash
mkdir data
cd data/

#download data
#The wget option -O specifies a file to which the documents is written, and here we use -, meaning it will written to standard output and piped to tar and the tar #flag -x enables extraction of archive files and -z decompresses, compressed archive files created by gzip

wget -c https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O - | tar -xz
cd ..
```
```bash
mkdir output
cd output
mkdir images 
mkdir timings
cd ..
```
```bash
mkdir bin
cd bin/
```
```bash
conda activate seurat
rstudio
```

open file: R.packages.R, install the packages \
open seurat_tutorial.Rmd
