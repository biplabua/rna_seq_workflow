# Snakemake workflow for differential gene expression between two GEO data sets of your choice

In order to run this pipeline install bioconda by following the instruction in the link below.

https://bioconda.github.io/user/install.html

Install snakemake by following the instruction in the link below.
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-pip

After installing snakemake run following command to run the pipeline. 

```conda activate snakemake```

```snakemake --core 1```

# To analyze geo data of your interest
Remove the default link and add link for the transcript sequence in the ```congif.yaml``` file.

If the data is pair end set ```paired=True``` in the ```congif.yaml``` file, otherwise set ```paired=False```. 

Remove the default sra id in ```meta-data.csv``` file and add your own sra id number.

Then run

```snakemake --core 1```
