# Fallstudie



## Get started


### Conda
**Create an environment from YAML file:**
 
copy seutat.yml file in project folder \
inside this folder:

```bash
conda env create --file seurat.yml
```

```bash
mkdir bin
```

```bash
mv <path to seurat_tutorial.Rmd>/seurat_tutorial.Rmd ./bin/
```

```bash
mkdir data
cd data

#download data
#The wget option -O specifies a file to which the documents is written, and here we use -, meaning it will written to standard output and piped to tar and the tar #flag -x enables extraction of archive files and -z decompresses, compressed archive files created by gzip

wget -c https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O - | tar -xz
##for our new dataset, do we need the raw or the filtered??
wget -c https://cf.10xgenomics.com/samples/cell exp/4.0.0/Targeted_SC3v3_Human_Glioblastoma_Neuroscience/Targeted_SC3v3_Human_Glioblastoma_Neuroscience_raw_feature_bc_matrix.tar.gz -O - | tar -xz
wget -c https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Targeted_SC3v3_Human_Glioblastoma_Neuroscience/Targeted_SC3v3_Human_Glioblastoma_Neuroscience_filtered_feature_bc_matrix.tar.gz -O - | tar -xz
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
conda activate seurat
rstudio
```

open seurat_tutorial.Rmd

Maybe you run in troubles, it could be helpful to delete the environment and create again like in step1

```bash
conda env remove -n seurat
```
### Docker

inside the  project folder:


```bash
mkdir bin
```

```bash
mv <path to seurat_tutorial.Rmd>/seurat_tutorial.Rmd ./bin/
```

```bash
mkdir data
cd data

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
sudo docker run --rm -d -p 8787:8787 -e PASSWORD=password -v $PWD:/home/rstudio/ seurat
```

open in your webbrowser: \
http://localhost:8787/ 

Username : rstudio \
Password: password 


open seurat_tutorial.Rmd


*stop running container*

```bash
sudo docker ps
sudo docker stop <container name>
```
