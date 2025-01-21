## Introduction
This Python script is an adaptation of [single-cell gene expression depth script](https://github.com/genome/scrna_mutations/tree/master/gex-depth-position) developed as part of [Petti et al. (2019)](https://www.nature.com/articles/s41467-019-11591-1).

The code has been reimplemented in Python3 with changes to inputs (see below) and reports a different metric than the original script; instead of reporting the number of unique reads (CB:UMI) per position, the proportion of cells with >= 1 read at each position is reported (#CBs/position). 

Please note that this script is a work in progress (some warnings may get thrown when running the script) so please report any issues/suggestions.

## Inputs
Please specify the location of the following required inputs in the `scripts/script_parameters.tsv` file:

| Input      | Notes                                                                                                                                                                                                                                                    |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| bam_file   | from Cell Ranger output (please make sure BAM index is also present in the same directory)                                                                                                                                                               |
| barcodes   | list of 10x barcodes from Cell Ranger, tsv file with one barcode per line                                                                                                                                                                                |
| genes      | list of genes to plot (please use HGNC gene symbols e.g. GAPDH). Can also contain optional mutation data to plot in 'mutations' columns. See tutorial below for example. Only expressed mutations i.e. those located in exonic positions should be input |
| method     | please specify data type, either "10x" or "Kx" (PacBio Kinnex scRNA-seq)                                                                                                                                                                                 |
| sample     | name of sample (will be included in output file names)                                                                                                                                                                                                   |

## Outputs
The following files will be generated in the `output` directory after running the script:

| Output                           | Notes                                                                                      |
| -------------------------------- | -------------------------------------------------------------------------------------------|
| sample_name_coverage_plots.jpg   | coverage plots generated with matplotlib - see tutorial below for example                  |
| coverage_results.csv             | dataframe containing gene coverage data - see table below for explanation                  |
| ensembl_genes.csv                | results from ENSEMBL REST API query - needed to ascertain exonic regions of input genes    |

### coverage_results.csv
The output `output/coverage_results.csv` takes the following format:

| genomic_position                                                                       | proportial_coverage                                    | percent_gene_expression                                                                  | method                                         | sample         | gene         |
| -------------------------------------------------------------------------------------- | ------------------------------------------------------ | ---------------------------------------------------------------------------------------- | ---------------------------------------------- | -------------- | ------------ |
| genomic coordinates of exonic regions of the canonical transcript of the specific gene | proportion of cells at position with at least one read | % of cells with expression of gene (will be the same for all positions of the same gene) | which technology was defined in the input file | name of sample | name of gene |

## 10x data tutorial
Sample files are included in this GitHub repo to generate 2 coverage plots (GAPDH and NPM1) for an AML 10x scRNA-seq sample. Please note that the original files were subset to only include data for GAPDH and NPM1 so no other genes can be used as input. To run the tutorial, please follow the commands below:

### Step 1: Build and activate conda environment
`conda env create -f scripts/sc_gene_coverage.yml`

`conda activate sc_gene_coverage`

### Step 2: Run gene_coverage.py script
`python3 scripts/gene_coverage.py`

### Example Output
The tutorial above should result in the generation of the plot below in the `output` directory. The red dashed line in the NPM1 plot indicates the position of the mutation speicified in the `input/sample_genes.tsv` file (sample_name_coverage_plots.jpg)
![Coverage plots from sample data](https://github.com/modalaigh/Bern_PacBio_Data_Day_25/blob/main/sc_gene_coverage/tutorial_results/sample_name_coverage_plots.jpg)

## 10x vs Kinnex data tutorial 
If you are interested in comparing coverage for dual sequenced samples, please see the [dual_plot](https://github.com/modalaigh/Bern_PacBio_Data_Day_25/tree/main/ds_sc_gene_coverage) branch for a similar tutorial
