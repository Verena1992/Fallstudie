# Fallstudie

conda --version # conda 4.10.3


an environment with seurat v4 was created:

conda create --name seurat -c conda-forge -c bioconda r-seurat=4*

and the environment exported to a YAML file:

conda activate seurat
conda env export --name seurat  > seurat.yml
conda deactivate



conda create --name seurat r-base=4.1.0
conda activate seurat
conda install r-ggplot2=3.3.5 r-patchwork=1.1.1 r-seuratObject=4.0.2 r-seurat=4.0.4 r-dplyr=1.0.7
conda install r-rmarkdown r-knitr r-xfun r-htmltools




## Get started

### Create an environment from YAML file:
 
 copy seutat.yml file in project folder \
 inside this folder:

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

#Docker

##Installation in Linux

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
docker --version
```
Docker version 20.10.7, build 20.10.7-0ubuntu5~20.04.2




