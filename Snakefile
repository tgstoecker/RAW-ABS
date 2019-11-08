def trimmomatic_inputs(wildcards):
    if (seq_type == "pe"):
        return expand("{reads}_{strand}.fastq", strand=["R1", "R2"], reads=wildcards.reads)
    elif (seq_type == "se"):
        return expand("{reads}.fastq", reads=wildcards.reads)

#modify trimmomatic like this - but I have to call two different wrappers...
rule bowtie2:
    input:
        bowtie2_inputs,
        index=bowtie2_index
    output:
        sam="{reads}_bowtie2.sam"
    run:
        if seq_type == "pe":
            shell("bowtie2 -x {input.index} -1 {input.forward} -2 {input.reverse} -S {output.sam}")
        elif seq_type == "se":
            shell("bowtie2 -x {input.index} -U {input.reads} -S {output.sam}")
    
    
    
rule fastqc_RAW:
    input:
        "reads/{sample}.fastq"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "v0.41.0/bio/fastqc"
    
    
rule trimmomatic_PE:
    input:
        r1="reads/{sample}.1.fastq.gz",
        r2="reads/{sample}.2.fastq.gz"
    output:
        r1="trimmed/{sample}.1.fastq.gz",
        r2="trimmed/{sample}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="trimmed/{sample}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        32
    wrapper:
        "0.40.2/bio/trimmomatic/pe"


rule trimmomatic_SE:
    input:
        "reads/{sample}.fastq.gz"  # input and output can be uncompressed or compressed
    output:
        "trimmed/{sample}.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        # optional compression levels from -0 to -9 and -11
        compression_level="-9"
    threads:
        32
    wrapper:
        "0.40.2/bio/trimmomatic/se"



rule fastqc_TRIMMED:
    input:
        "reads/{sample}.fastq"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "v0.41.0/bio/fastqc"
        

rule star_index:
    input:
        fasta = "{genome}.fasta"
    output:
        directory("{genome}")
    message:
        "Testing STAR index"
    threads:
        1
    params:
        extra = ""
    log:
        "logs/star_index_{genome}.log"
    wrapper:
        "0.39.0/bio/star/index"


rule star_PE_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = ["reads/{sample}_R1.1.fastq", "reads/{sample}_R1.2.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = ["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"] #optional
    output:
        # see STAR manual for additional output files
        "star/pe/{sample}/Aligned.out.sam"
    log:
        "logs/star/pe/{sample}.log"
    params:
        # path to STAR reference genome index
        index="index",
        # optional parameters
        extra=""
    threads: 8
    wrapper:
        "0.39.0/bio/star/align"


rule star_SE:
    input:
        fq1 = "reads/{sample}_R1.1.fastq"
    output:
        # see STAR manual for additional output files
        "star/{sample}/Aligned.out.sam"
    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index="index",
        # optional parameters
        extra=""
    threads: 8
    wrapper:
        "0.39.0/bio/star/align"
