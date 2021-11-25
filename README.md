# Fallstudie



## Get started


### Conda
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

open seurat_tutorial.Rmd

### Docker

inside the  project folder:

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
sudo docker run --rm -d -p 8787:8787 -e PASSWORD=password -v $PWD:/home/rstudio/ seurat
```
http://localhost:8787/
Username : rstudio
Password: password


open seurat_tutorial.Rmd



*stop running container*

```bash
sudo docker ps
sudo docker stop <container name>
```











