##### load config and define samples ####

SAMPLE,  = glob_wildcards("rawreads/{sample}.fq.gz")
SAMPLES, = glob_wildcards("rawreads/{sample}_1.fq.gz")
configfile: "config.yaml"

#### target rules ####

rule all:
    input:
#        expand("fastqc/raw/{sample}_{paired}_fastqc.html", sample=SAMPLES, paired=[1, 2]),
#        expand("fastqc/raw/{sample}_{paired}_fastqc.zip", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.1.fq.gz", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.2.fq.gz", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.1.unpaired.fq.gz", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.2.unpaired.fq.gz", sample=SAMPLES, paired=[1, 2]),
##        expand("star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
#        expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}.dedup.txt", sample=SAMPLES),
        expand("star/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand("removed_duplicates_alignments/{sample}.dedup.bam.bai", sample=SAMPLES),
        "multiqc/multiqc.html",
         directory(expand("FGS/{genome}", genome=config["genome"]))


rule STAR_index:
    input:
        fasta = expand("FGS/{genome}.fa", genome=config["genome"])
    output:
        directory(expand("FGS/{genome}", genome=config["genome"]))
    message:
        "Creating STAR index"
    params:
        extra = "",
        gtf = expand("FGS/{annotation}.gtf", annotation=config["annotation"]),
        threads= config["threads_star_index"],
        length = config["read_length_star_index"]
    log:
        expand("logs/star_index_{genome}.log", genome=config["genome"])
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {params.threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fasta} '
            '--sjdbGTFfile {params.gtf} '
            '--sjdbOverhang {params.length}'


##########################################################################################
##########################################################################################
## For Single End Reads -SE
if config["sequencing_type"] == "single_end":
    rule fastqc_SE:
        input:
            expand("rawreads/{sample}.fq.gz", sample=SAMPLE)
        output:
            html="fastqc/raw/{sample}_fastqc.html",
            zip="fastqc/raw/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        #params: "-t 2 --quiet"
        log:
            "logs/fastqc/raw/{sample}.log"
        shell:
            "fastqc -t 2 --quiet {input} -o fastqc/raw/"


    rule trimmomatic_SE:
        input:
            "rawreads/{sample}.fq.gz"
        output:
            "trimmed/{sample}.fq.gz"
        log:
            "logs/trimmomatic/{sample}.log"
        params:
            #  list of trimmers (see manual)
            #trimmer=["SLIDINGWINDOW:4:15"],
            trimmer={config["trim.options"]},
            # optional parameters
            extra="",
            compression_level="-9"
        threads: 4
        shell:
            "trimmomatic SE -threads {threads} {input} {output} {params.trimmer}"


    rule trimmed_fastqc_SE:
        input:
            expand("trimmed/{sample}.fq.gz", sample=SAMPLE)
        output:
            html="fastqc/trimmed/{sample}_fastqc.html",
            zip="fastqc/trimmed/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        params: "-t 2 --quiet"
        log:
            "logs/fastqc/trimmed/{sample}.log"
        shell:
            "fastqc {params} {input} -o fastqc/trimmed/"


    rule STAR_SE:
        input:
            file = "trimmed/{sample}.fq.gz",
            dir = expand("FGS/{genome}", genome=config["genome"])
        output:
            # see STAR manual for additional output files - 
            "star/{sample}.Aligned.sortedByCoord.out.bam"
        params:
            name = "star/{sample}.",
            threads = config["threads_star"]
        log:
            "logs/star/{sample}.log"
        shell:
            'STAR --runThreadN {params.threads} '
            '--genomeDir {input.dir} '
            '--readFilesIn {input.file} '
            '--readFilesCommand zcat '
            '--outSAMtype BAM SortedByCoordinate '
            '--outFileNamePrefix {params.name} '
            '--outBAMsortingThreadN {params.threads}'


    rule multiqc_SE:
        input:
            expand("fastqc/raw/{sample}_fastqc.zip", sample=SAMPLE),
            expand("logs/trimmomatic/{sample}.log", sample=SAMPLE),
            expand("fastqc/trimmed/{sample}_fastqc.zip", sample=SAMPLE),
            expand("star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLE),
            expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLE)
        output:
            "multiqc/multiqc.html"
        params:
            "-ip"
        log:
            "logs/multiqc.log"
        wrapper:
            "0.42.0/bio/multiqc"


##########################################################################################
##########################################################################################
### For Paired End Reads -PE

if config["sequencing_type"] == "paired_end":
    rule fastqc_PE:
        input:
            expand("rawreads/{sample}_{paired}.fq.gz", sample=SAMPLES, paired=[1, 2])
        output:
            html="fastqc/raw/{sample}_{paired}_fastqc.html",
            zip="fastqc/raw/{sample}_{paired}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        params: "-t 2 --quiet"
        log:
            "logs/fastqc/raw/{sample}_{paired}.log"
        shell:
            "fastqc {params} {input} -o fastqc/raw/"


    rule trimmomatic_PE:
        input:
            r1="rawreads/{sample}_1.fq.gz",
            r2="rawreads/{sample}_2.fq.gz"
        output:
            r1="trimmed/{sample}.forward_paired.fq.gz",
            r2="trimmed/{sample}.reverse_paired.fq.gz",
            # reads where trimming entirely removed one of the mates
            r1_unpaired="trimmed/{sample}.forward_unpaired.fq.gz",
            r2_unpaired="trimmed/{sample}.reverse_unpaired.fq.gz"
        log:
            "logs/trimmomatic/{sample}.log"
        params:
            #  list of trimmers (see manual)
            #trimmer=["SLIDINGWINDOW:4:15"],
            trimmer={config["trim.options"]},
            # optional parameters
            extra="",
            compression_level="-9"
        threads: 4
        wrapper:
            "0.42.0/bio/trimmomatic/pe"


    rule trimmed_fastqc_PE:
        input:
            expand("trimmed/{sample}.{mate}_paired.fq.gz", sample=SAMPLES, mate=["forward", "reverse"])
        output:
            html="fastqc/trimmed/{sample}.{mate}_paired_fastqc.html",
            zip="fastqc/trimmed/{sample}.{mate}_paired_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        params: "-t 2 --quiet"
        log:
            "logs/fastqc/trimmed/{sample}_{mate}.log"
        shell:
            "fastqc {params} {input} -o fastqc/trimmed/"



    rule STAR_PE:
        input:
            fq1 = "trimmed/{sample}.forward_paired.fq.gz",
            fq2 = "trimmed/{sample}.reverse_paired.fq.gz",
            dir = expand("FGS/{genome}", genome=config["genome"])
        output:
            # see STAR manual for additional output files -
            "star/{sample}.Aligned.sortedByCoord.out.bam"
        params:
            name = "star/{sample}.",
            threads = config["threads_star"]
        log:
            "logs/star/{sample}.log"
        params:
            # path to STAR reference genome index
            #index=lambda wildcards: expand("FGS/{genome}", genome=config["genome"]),
            #threads = config["threads_star"]
        shell:
            'STAR --runThreadN {params.threads} '
            '--genomeDir {input.dir} '
            '--readFilesIn {input.fq1} {input.fq2} '
            '--readFilesCommand zcat '
            '--outSAMtype BAM SortedByCoordinate '
            '--outFileNamePrefix {params.name} '
            '--outBAMsortingThreadN {params.threads}'


    rule multiqc_PE:
        input:
            expand("fastqc/raw/{sample}_{paired}_fastqc.zip", sample=SAMPLES, paired=[1, 2]),
            expand("logs/trimmomatic/{sample}.log", sample=SAMPLES),
            expand("fastqc/trimmed/{sample}.{mate}_paired_fastqc.zip", sample=SAMPLES, mate=["forward", "reverse"]),
            expand("star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
            expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES)
        output:
            "multiqc/multiqc.html"
        params:
            "-ip"
        log:
            "logs/multiqc.log"
        wrapper:
            "0.42.0/bio/multiqc"


##########################################################################################
##########################################################################################
### Continuation regardless of sequencing type


rule index_sorted_bams_with_dups:
    input:
        "star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "star/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params:
        threads = config["threads_index_sorted_bams_with_dups"]
    shell:
        "samtools index -@ {params.threads} {input}"


rule remove_duplicates_picard:
    input:
         "star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
         bam="removed_duplicates_alignments/{sample}.dedup.bam",
         txt="removed_duplicates_alignments/{sample}.dedup.txt"
    log:
         "logs/picard/{sample}.dedup.log"
    shell:
         "picard MarkDuplicates I={input} O={output.bam} M={output.txt} REMOVE_DUPLICATES=true > {log} 2>&1"


rule index_sorted_bams_without_dups:
    input:
        "removed_duplicates_alignments/{sample}.dedup.bam"
    output:
        "removed_duplicates_alignments/{sample}.dedup.bam.bai"
    params:
        threads = config["threads_index_sorted_bams_without_dups"]
    shell:
        "samtools index -@ {params.threads} {input}"


rule featureCounts:
#The "-O" option is not necessary if your aligner is not junction-aware, 
#namely they don't report the exon-exon junctions in the CIGAR strings in your mapping results. 
#However, if your aligner is junction-aware, for example, HISAT or STAR, 
#then you HAVE to use the "-O" option, or you will lose all the reads or read-pairs that overlap with multiple exons.
#The "-O" option deals with the reads or read-pairs that overlaps with multiple exons or genes. 
#Say, a read is mapped to only one location, but there are 3 exons that all overlap with this location, 
#then you have to use "-O" to have this read counted (it contributes one count to each of the 3 exons). 
#If you don't use "-O", this read will be assigned to no exon because of the ambiguity.


#need to perform four different kinds
#with/without duplicates X gene/transcript level

## gene level
#featureCounts -T 8 -O -t exon -g gene_id -a $annotation/zea_mays.protein_coding.gtf -o $count/gene-level/total_file.count \
#B73_con_1_trimmed_sorted.bam \
#B73_con_2_trimmed_sorted.bam \
#B73_con_3_trimmed_sorted.bam \
#B73_con_4_trimmed_sorted.bam


    
