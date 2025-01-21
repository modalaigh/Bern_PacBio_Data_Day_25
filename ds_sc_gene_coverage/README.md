## Introduction
This is a very similar tutorial to the main tutorial but which focuses on comparing coverage of genes in cells which were sequenced with 2 technologies (e.g. 10x vs Kinnex). The genes/sample here are the same as the main tutorial (NPM1 and GAPDH in an AML relapse sample)

## Inputs
Please specify the location of the following required inputs in the `scripts/script_parameters.tsv` file. Note: some rows should contain csv values for your two data types - see notes below and example `scripts/script_parameters.tsv` file.

| Input      | Notes                                                                                                                                                                                                                                                    |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| bam_file*  | from Cell Ranger/SMRT Link output (please make sure BAM index is also present in the same directory)                                                                                                                                                     |
| barcodes*  | list of 10x/Kinnex barcodes from Cell Ranger/SMRT Link, tsv file with one barcode per line                                                                                                                                                               |
| genes      | list of genes to plot (please use HGNC gene symbols e.g. GAPDH). Can also contain optional mutation data to plot in 'mutations' columns. See tutorial below for example. Only expressed mutations i.e. those located in exonic positions should be input |
| method     | please specify order of data types, either "10x,Kx" or "Kx,10x". Note: other files (those with *) should be specified in the same order                                                                                                                  |
| sample     | name of sample (will be included in output file names)                                                                                                                                                                                                   |

> [!NOTE]  
> The script identifies cells which were sequenced with both technologies (dual-sequenced cells) and looks at the gene coverage of only these cells. That is why there are only 4456 cells visualised in the graph compared with 7051 in the main 10x-only tutorial.

## Outputs
The following files will be generated in the `output` directory after running the script:

| Output                           | Notes                                                                                      |
| -------------------------------- | -------------------------------------------------------------------------------------------|
| sample_name_coverage_plots.jpg   | coverage plots generated with matplotlib - see tutorial below for example                  |
| coverage_results.csv             | dataframe containing gene coverage data - see table below for explanation                  |
| ensembl_genes.csv                | results from ENSEMBL REST API query - needed to ascertain exonic regions of input genes    |

### coverage_results.csv
The output `output\coverage_results.csv` takes the following format:

| genomic_position                                                                       | proportial_coverage                                    | percent_gene_expression                                                                                             | method                            | sample         | gene         |
| -------------------------------------------------------------------------------------- | ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------- | --------------------------------- | -------------- | ------------ |
| genomic coordinates of exonic regions of the canonical transcript of the specific gene | proportion of cells at position with at least one read | % of cells with expression of gene (will be the same for all positions of the same gene for each sequecning method) | the sequencing method (10x vs Kx) | name of sample | name of gene |

## Tutorial
Sample files are included in this GitHub repo to generate 2 coverage plots (GAPDH and NPM1) for a dual-sequenced AML scRNA-seq sample (10x and Kinnex). Please note that the original files were subset to only include data for GAPDH and NPM1 so no other genes can be used as input. To run the tutorial, please follow the commands below:

### Step 1: Build and activate conda environment
`conda env create -f scripts/sc_gene_coverage.yml`

`conda activate sc_gene_coverage`

### Step 2: Run dual_gene_coverage.py script
`python3 script/gene_coverage.py`

### Example output
The tutorial above should result in the generation of the plot below. The red dashed line in the NPM1 plot indicates the position of the mutation speicified in the `script/sample_genes.tsv` file. The other colours correspond with the sequencing method (see legend)

![Coverage plots from sample data](https://github.com/modalaigh/Bern_PacBio_Data_Day_25/blob/main/ds_sc_gene_coverage/tutorial_results/sample_name_coverage_plots.jpg)
