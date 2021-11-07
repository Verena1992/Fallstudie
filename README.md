# Fallstudie

an environment with seurat v4 was created:

conda create --name seurat -c conda-forge -c bioconda r-seurat=4*

and the environment exported to a YAML file:

conda activate seurat
conda env export --name seurat  > seurat.yml
conda deactivate



Create an environment from YAML file:

conda env create --file seurat.yml



mkdir data
cd data/
#The wget option -O specifies a file to which the documents is written, and here we use -, meaning it will written to standard output and piped to tar and the tar #flag -x enables extraction of archive files and -z decompresses, compressed archive files created by gzip

wget -c https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O - | tar -xz
cd ..

#make data read only
chmod -wx data/*
mkdir output
cd output
mkdir images 
mkdir timings

cd ..


mkdir bin
cd bin/
conda activate seurat
rstudio
(open seurat_tutorial.rmd)
