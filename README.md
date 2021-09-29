


# Requirements

- python 3.6.8
- Cell Ranger v5.0.1

## Python packages
- pepars 
- scrapi
- scvi (custom fork: https://github.com/dibidave/scvi-tools)
- pandas
- numpy
- skbio
- statsmodels
- scipy

# Installing

To install the included package, you must clone it and install it locally:
```
git clone git@github.com:GradinaruLab/aavomics.git
cd aavomics
pip install --upgrade .
```

# Setup

Data for this work can be reproduced in two different ways:

1. From preprocessed gene count matrices (this only requires about 10GB of space)
2. From raw FASTQ files (note - this means you will need heavy compute time to align the samples, and about 2TB of storage for the reads)

Regardless of which route you take, however, you must do the following:

### 1. Set your data path environment variable

This package and all scripts rely on a fixed structure of files. To make sure this is consistent across scripts, we provide an environment variable ```AAVOMICS_PATH```, which, when set, all calls to the aavomics package will use to look for data. By default, this path will be a directory "data" relative to wherever you are running - therefore, to use the same data across scripts, you must set the environment variable. You can do this through bash, or you can do it in python (before importing aavomics) by:
```
os.environ["AAVOMICS_PATH"] = MY_PATH
```
### 2. Download metadata and processed data

The aavomics package relies on the sample and other metadata to be stored in CSV files in a custom format. This data maps the samples to their associated reads, metadata, read sets, and more. This is necessary to run any of the figure generation scripts. To download the metadata, run the "Download the database files" section of ```examples/Download Data.ipynb```.

## From processed data

If you want to skip preprocessing and go straight to the final gene counts, you can skip right to the figure generation scripts after downloading the metadata and processed data.

## From raw data

### 1. Download raw data

Raw FASTQ files can be downloaded from SRA project https://www.ncbi.nlm.nih.gov/bioproject/758711, and must be placed within their respective cell_set folder. Transcriptome reads go in ```AAVOMICS_PATH/cell_sets/{CELL_SET_NAME}/transcriptome/reads``` and amplified viral transcript reads go in ```AAVOMICS_PATH/cell_sets/{CELL_SET_NAME}/virus/reads```

### 2. Generate custom references

To align the reads against a reference that includes the AAV gene, we create a custom reference. An example of this is in: ```examples/alignment/Create Custom Reference.ipynb```

### 3. Align samples

We perform multiple alignments - one to match the Allen brain atlas reference, one with a custom AAV reference, and one with a standard reference for cell typing (to avoid virus counts from being included in cell typing). Alignments are run as slurm jobs that run Cell Ranger. An example of how to run alignment is in ```examples/alignment/Align Transcriptomes.ipynb```, which is a Python script that creates modified versions of the slurm job script in ```examples/bin/cellranger_count.sh```.

### 4. Convert to h5ad

We convert all Cell Ranger h5 files to anndata files for convenience. An example of doing this is in ```examples/alignment/Convert Cell Ranger H5 to Anndata.ipynb```

### 5. Align virus reads

We align virus reads using a custom wrapper around the Striped Smith–Waterman alignment from ```skbio.alignment```. Since we expect perfect matches with the exception of possible barcode insertions, we align all reads against our expected plasmid templates, and then save only the alignment CIGAR string and any associated insertions. An example of this is in ```examples/alignment/Align Virus Reads.ipynb```

### 6. Demultiplex viral transcripts

After aligning viral reads, we use these amplified reads to demultiplex reads that align to the custom AAV gene from our transcriptome. For each aligned read in the transcriptome, the variant gets called if that variant barcode or template is present for that cell barcode UMI higher than any other viral variant. This is done in ```examples/alignment/Demultiplex Viral Counts.ipynb```

### 7. Cell Typing

To reproduce our cell typing pipeline, you can start by running the ```examples/cell_typing/Droplet Classification.ipynb``` notebook, which trains a model on a subsampling of the cells across all cell sets, and predicts whether each droplet is empty, a multiplet, a neuron or a non-neuron, based on whether the cluster expresses the canonical marker genes.

After droplets have been classified, you can type the neurons by running ```examples/cell_typing/Neuron Subtyping Allen Reference.ipynb``` (note: you must first download and filter the Allen dataset, i.e. via ```examples/cell_typing/Filter Allen Dataset.ipynb```.

You can then type non-neurons by running ```examples/cell_typing/Non-Neuron Subtyping.ipynb```.

Finally, to combine the cell types, you can run ```examples/Combine Neuron and Non-Neuron Subtyping.ipynb``` and then to save them out to a single anndata for analysis you can run ```examples/cell_typing/Make Combined Anndata.ipynb```.

### 8. Analysis

You can also reproduce our analysis scripts, which convert the viral counts into estimates of viral transduction by cell type or cluster, by running ```examples/analysis/Estimate Viral Transduction by Cell Type.ipynb```.

Data is prepared in a format that is amenable for DESeq2 input in ```examples/figures/Figure 5B,C.ipynb``` and the corresponding .csv files files are used as in input to ```examples/figures/Figure 5B.R```. Files generated by DESeq2 are stored in the "deseq2_out" folder and the rest of the ```examples/figures/Figure 5B,C.ipynb``` script is used to generate diffrential expression heatmaps across 3DPI and 25DPI conditions for each cell type, and the histogram plot indicating number of differentially expressed genes for each cell type. 
Script used to compare gene expression of gene of interest between 3DPI and 25DPI for cell of interest is in ```examples/figures/Figure 5E.ipynb```




# Figure reproduction

All figures in the paper can be produced by running the following scripts off of the processed data:

- Figure 2C,D: ```examples/figures/PHP.eB CAP-B10 Coinjection IHC vs scRNA-seq.ipynb```
- Figure 3: ```examples/figures/CAP-B10 vs PHP.eB Neuron Tropism.ipynb```
- Figure 4B: ```examples/figures/Non-Neuron t-SNE.ipynb```
- Figure 4C: ```examples/figures/Non Neuron Marker Gene Dot Plot.ipynb```
- Figure 4D: ```examples/figures/Barcode vs Cargo Variability.ipynb```
- Figure 4E: ```examples/figures/PHP-V1 vs PHP-eB Major Cell Type Barcode vs Cargo.ipynb```
- Figure 4F: ```examples/figures/PHP-eB vs PHP-V1 Delta Fraction Transduced.ipynb.ipynb```
- Figure 5.A: ```Figure 5A.ipynb```
- Figure 5.B: ```Figure_5B.R```
- Figure 5.B: ```Figure 5B,C.ipynb```
- Figure 5.C: ```Figure 5B,C.ipynb```
- Figure 5.E: ```Figure 5E.ipynb```
- Figure 6: ```examples/figures/BC4 Tropism Dot Plot.ipynb```
- Supplemental Figure 2A: ```examples/figures/Threshold Transduction Rates.ipynb```
- Supplemental Figure 2B: ```examples/Mean Transcripts vs Cumulative Probability.ipynb```
- Supplemental Figure 2D: ```examples/analysis/Counting Simulation.ipynb``` (note: requires processing from raw data to get access to debris distributions)
- Supplemental Figure 2E: ```examples/figures/Transduction Rate vs Transcripts per Cell.ipynb```
- Supplemental Figure 3A,B: ```examples/figures/Cell Ranger Debris Example.ipynb```
- Supplemental Figure 3C: ```examples/figures/Missing Multiplets Example.ipynb```
- Supplemental Figure 4B: ```examples/figures/Cell Type Sunburst.iypnb```
- Supplemental Figure 4C: ```examples/figures/Cell Type Box Plots.ipynb```
- Supplemental Figure 4D: ```examples/figures/Empty Droplet QC.ipynb```
- Supplemental Figure 5A: ```examples/figures/BC Expression Variability.ipynb```
- Supplemental Figure 5B: ```examples/figures/Fluorophore Variability.ipynb```
- Supplemental Figure 5C: ```examples/figures/Fluorophore Variability Vascular Tropism.ipynb```
- Supplemental Figure 6A: ```examples/figures/Intra vs Inter Sample Variance.ipynb```
- Supplemental Figure 6B: ```examples/figures/Percent of Transduced vs Percent Total.ipynb```
- Supplemental Figure 7: ```examples/figures/Subtype Transduction Rate Comparisons.ipynb```

