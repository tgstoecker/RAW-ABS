rule fastqc:
    input:
        "rawreads/{sample}.{paired}.fastq"
    output:
        html="qc/fastqc/{sample}.{paired}.html",
        zip="qc/fastqc/{sample}.{paired}._fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.{paired}.log"
    wrapper:
        "v0.41.0/bio/fastqc"
