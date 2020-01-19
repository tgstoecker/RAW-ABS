# RAW-ABS - RNAseq Analysis Workflow, Alignment BaSed 
Alignment based workflow for Illumina Short reads


alignment based (could also have second version that functions "alignment free" so e.g. salmon - however since these are all based on aligning to a transcriptome index -> sth. like my trinity workflow or stringtie would be necessary excluding cases, where a reference actually exist and we choose not to use it for some reason..)

Install the Python 3 version of Miniconda.
you can get it here: https://docs.conda.io/en/latest/miniconda.html

Answer yes to the question whether conda shall be put into your PATH.
For detailed options concerning conda/bioconda see:

Then, you can install Snakemake with

`conda install -c bioconda -c conda-forge snakemake`

Preparing a working directory
First, create a new directory and change into that directory in your terminal.

Download/Clone the current release of the MuWU pipeline into the directory.

The included environment.yaml file can be used to install all required software into an isolated Conda environment with a name of your choice - in the following we will call it "snakemake-MuWU":

`conda env create --name RAW-ABS --file environment.yaml`

Activating the environment
To activate the snakemake-tutorial environment, execute

`conda activate RAW-ABS`

Now you can use the installed tools and our workflow without any software dependency issues.
For detailed options of snakemake see: 

Should you want to remove the conda environment, execute
`conda env remove -n RAW-ABS`


using Snakemake and an elaborate YAML or JSON file to control all options

just linking to the input files is enough - no need to copy!

sample data - phytophtara infestans
`wget ftp://ftp.ensemblgenomes.org/pub/protists/release-45/fasta/phytophthora_infestans/dna/Phytophthora_infestans.ASM14294v1.dna.toplevel.fa.gz`

`wget ftp://ftp.ensemblgenomes.org/pub/protists/release-45/gtf/phytophthora_infestans/Phytophthora_infestans.ASM14294v1.45.gtf.gz`

`fastq-dump --split-files https://sra-download.ncbi.nlm.nih.gov/traces/dra2/DRR/000156/DRR160421`

# What it does:
- FastQC on RawReads
- adapter removal (probably easiest use when I just include a directory where a adapter fasta should be put) /
- Trimming (switch to sth. better than trimmomatic?)
- FastQC on Trimmed Reads
- keep fastq_screen "database" and perform screening of all samples
- get rRNA contamination statistics using bbduk.sh (further stuff?)
- alignment/mapping: STAR (including index) - output directly to sorted bam with option to have sam/use sambamba?
- featureCounts on alignment files
- multiqc (remember e.g. gastq_screen multiqc output sucks!)
- (also complete R diff. exp. analysis?)
