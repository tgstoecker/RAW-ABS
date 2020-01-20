# RAW-ABS - RNAseq Analysis Workflow, Alignment BaSed 
# Alignment based workflow for Illumina Short reads


# Setup:
Install the Python 3 version of Miniconda.
you can get it here: https://docs.conda.io/en/latest/miniconda.html

Answer yes to the question whether conda shall be put into your PATH.
For detailed options concerning conda/bioconda see:

Then, you can install Snakemake with

`conda install -c bioconda -c conda-forge snakemake`

Preparing a working directory
First, create a new directory and change into that directory in your terminal.

Download/Clone the current release of the MuWU pipeline into the directory.

The included environment.yaml file can be used to install all required software into an isolated Conda environment with a name of your choice - in the following we will call it "RAW-ABS":

`conda env create --name RAW-ABS --file environment.yaml`

Activating the environment
To activate the snakemake-tutorial environment, execute

`conda activate RAW-ABS`

Now you can use the installed tools and our workflow without any software dependency issues.
For detailed options of snakemake see: 

Should you want to remove the conda environment, execute
`conda env remove -n RAW-ABS`

# Usage:
1) Move, copy or link fasta and gtf of your species into the FGS directory
2) Move, copy or link your gzipped fastq files to the rawreads directory
3) All options of the workflow can easily be controlled via the config.yaml file
  - rename your fastq files to follow the naming scheme: xxxx_1.fq.gz for PE reads!
4) When all this is done execute `snakemake -np` to check if the workflow works and perform a dry-run
5) TO start the workflow execute `snakemake --cores xx` and set the total amount of threads to be used

# What it does:
- Removal of rRNA reads via bbduk
- FastQC on rRNA depleted RawReads
- Trimming and FastQC on these trimmed Reads
- alignment/mapping via STAR (including index) - outputs directly to sorted bam
- featureCounts on alignment files
- multiqc
- all indexes for visualization software
