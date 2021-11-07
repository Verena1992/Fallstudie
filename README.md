# Fallstudie

an environment with seurat v4 was created:

conda create --name seurat -c conda-forge -c bioconda r-seurat=4*

and the environment exported to a YAML file:

conda activate seurat
conda env export --name seurat  > seurat.yml
conda deactivate



Create an environment from YAML file:

conda env create --file seurat.yml
