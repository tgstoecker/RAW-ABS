# RNA-seq-analysis-workflow-alignment-based
RNA-seq analysis workflow (alignment based)

alignment based (could also have second version that functions "alignment free" so e.g. salmon - however since these are all based on aligning to a transcriptome index -> sth. like my trinity workflow or stringtie would be necessary excluding cases, where a reference actually exist and we choose not to use it for some reason..)

using Snakemake and an elaborate YAML or JSON file to control all options

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
