##### load config and define samples ####

configfile: "config.yaml"
SAMPLES, = glob_wildcards(config['sample_names'])

#### target rules ####

rule all:
    input:
#        expand("fastqc/raw/{sample}_{paired}_fastqc.html", sample=SAMPLES, paired=[1, 2]),
#        expand("fastqc/raw/{sample}_{paired}_fastqc.zip", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.1.fq.gz", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.2.fq.gz", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.1.unpaired.fq.gz", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed/{sample}.2.unpaired.fq.gz", sample=SAMPLES, paired=[1, 2]),
#        expand("star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}.dedup.txt", sample=SAMPLES),
#        directory(expand("FGS/{genome}", genome=config["genome"])),
        "FGS/rRNA.gtf",
        "FGS/rRNA.fa",
        expand("star/{sample}.Aligned.sortedByCoord.out.bam.csi", sample=SAMPLES),
        expand("removed_duplicates_alignments/{sample}.dedup.bam.csi", sample=SAMPLES),
        "FC/with_dups_gene_level/counts.txt",
        "FC/with_dups_exon_level/counts.txt",
        "FC/without_dups_gene_level/counts.txt",
        "FC/without_dups_exon_level/counts.txt",
        "multiqc/multiqc.html"


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
        length = config["read_length_star_index"],
        size = config["limitGenomeGenerateRAM"]
    log:
        expand("logs/star_index_{genome}.log", genome=config["genome"])
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {params.threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--limitGenomeGenerateRAM {params.size}'
            '--genomeFastaFiles {input.fasta} '
            '--sjdbGTFfile {params.gtf} '
            '--sjdbOverhang {params.length}'


rule create_rRNA_gtf:
    input:
        gtf = expand("FGS/{annotation}.gtf", annotation=config["annotation"])
    output:
        "FGS/rRNA.gtf"
    shell:
        "grep 'rRNA' {input.gtf} > {output}"


rule create_rRNA_fasta:
    input:
        gtf = "FGS/rRNA.gtf",
        fa = expand("FGS/{genome}.fa", genome=config["genome"])
    output:
        fa = "FGS/rRNA.fa"
    shell:
        "bedtools getfasta -fi {input.fa} -bed {input.gtf} -fo {output.fa}"


##########################################################################################
##########################################################################################
### For Single End Reads -SE

if config["sequencing_type"] == "single_end":
    rule bbduk:
        input:
            reads = "rawreads/{sample}.fq.gz",
            fa = "FGS/rRNA.fa"
        output: 
            reads = "rawreads_rRNA_removed/{sample}.fq.gz",
            stats = "rawreads_rRNA_removed/{sample}.stats.txt"
        params:
            threads = config["threads_bbduk"],
            kmer_length = config["kmer_length"] 
        shell: 
            "bbduk.sh k={params.kmer_length} threads={params.threads} in={input.reads} out={output.reads} ref={input.fa} stats={output.stats}"


    rule fastqc:
        input:
            "rawreads_rRNA_removed/{sample}.fq.gz"
        output:
            html="fastqc/raw/{sample}_fastqc.html",
            zip="fastqc/raw/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        #params: "-t 2 --quiet"
        log:
            "logs/fastqc/raw/{sample}.log"
        shell:
            "fastqc -t 2 --quiet {input} -o fastqc/raw/"


    rule trimmomatic:
        input:
            "rawreads_rRNA_removed/{sample}.fq.gz"
        output:
            "trimmed/{sample}.fq.gz"
        log:
            "logs/trimmomatic/{sample}.log"
        params:
            trimmer={config["trim.options"]},
            # optional parameters
            extra="",
            compression_level="-9"
        threads: 4
        shell:
            "trimmomatic SE -threads {threads} {input} {output} {params.trimmer}"


    rule trimmed_fastqc:
        input:
            "trimmed/{sample}.fq.gz"
        output:
            html="fastqc/trimmed/{sample}_fastqc.html",
            zip="fastqc/trimmed/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        params: "-t 2 --quiet"
        log:
            "logs/fastqc/trimmed/{sample}.log"
        shell:
            "fastqc {params} {input} -o fastqc/trimmed/"


    rule STAR:
        input:
            file = "trimmed/{sample}.fq.gz",
            dir = expand("FGS/{genome}", genome=config["genome"])
        output:
            # see STAR manual for additional output files - 
            "star/{sample}.Aligned.sortedByCoord.out.bam"
        params:
            name = "star/{sample}.",
            threads = config["threads_star"],
            outReadsUnmapped = config["outReadsUnmapped"],
            outSAMattributes = config["outSAMattributes"],
            outSJfilterOverhangMin = config["outSJfilterOverhangMin"],
            outFilterMultimapNmax = config["outFilterMultimapNmax"],
            outFilterScoreMin = config["outFilterScoreMin"],
            outFilterMatchNminOverLread = config["outFilterMatchNminOverLread"],
            outFilterMismatchNmax = config["outFilterMismatchNmax"],
            outFilterMismatchNoverLmax = config["outFilterMismatchNoverLmax"],
            alignIntronMin = config["alignIntronMin"],
            alignIntronMax = config["alignIntronMax"],
            alignMatesGapMax = config["alignMatesGapMax"],
            alignSJoverhangMin = config["alignSJoverhangMin"],
            alignSJDBoverhangMin = config["alignSJDBoverhangMin"],
            alignSoftClipAtReferenceEnds = config["alignSoftClipAtReferenceEnds"],
            chimSegmentMin = config["chimSegmentMin"],
            chimScoreMin = config["chimScoreMin"],
            chimScoreSeparation = config["chimScoreSeparation"],
            chimJunctionOverhangMin = config["chimJunctionOverhangMin"],
            quantMode = config["quantMode"],
            twopassMode = config["twopassMode"],
            chimOutType = config["chimOutType"],
            outSAMtype = config["outSAMtype"],
            sjdbOverhang = config["sjdbOverhang"],
            outFilterMultimapScoreRange = config["outFilterMultimapScoreRange"]
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
            '--outReadsUnmapped {params.outReadsUnmapped}'
            '--outSAMattributes {params.outSAMattributes}'
            '--outSJfilterOverhangMin {params.outSJfilterOverhangMin}'
            '--outFilterMultimapNmax {params.outFilterMultimapNmax}'
            '--outFilterScoreMin {params.outFilterScoreMin}'
            '--outFilterMatchNminOverLread {params.outFilterMatchNminOverLread}'
            '--outFilterMismatchNmax {params.outFilterMismatchNmax}'
            '--outFilterMismatchNoverLmax {params.outFilterMismatchNoverLmax}'
            '--alignIntronMin {params.alignIntronMin}'
            '--alignIntronMax {params.alignIntronMax}'
            '--alignMatesGapMax {params.alignMatesGapMax}'
            '--alignSJoverhangMin {params.alignSJoverhangMin}'
            '--alignSJDBoverhangMin {params.alignSJDBoverhangMin}'
            '--alignSoftClipAtReferenceEnds {params.alignSoftClipAtReferenceEnds}'
            '--chimSegmentMin {params.chimSegmentMin}'
            '--chimScoreMin {params.chimScoreMin}'
            '--chimScoreSeparation {params.chimScoreSeparation}'
            '--chimJunctionOverhangMin {params.chimJunctionOverhangMin}'
            '--quantMode {params.quantMode}'
            '--twopassMode {params.twopassMode}'
            '--chimOutType {params.chimOutType}'
            '--outSAMtype {params.outSAMtype}'
            '--sjdbOverhang {params.sjdbOverhang}'
            '--outFilterMultimapScoreRange {params.outFilterMultimapScoreRange}'


    rule multiqc:
        input:
            expand("rawreads_rRNA_removed/{sample}.stats.txt", sample=SAMPLES),
            expand("fastqc/raw/{sample}_fastqc.zip", sample=SAMPLES),
            expand("logs/trimmomatic/{sample}.log", sample=SAMPLES),
            expand("fastqc/trimmed/{sample}_fastqc.zip", sample=SAMPLES),
            expand("star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
            expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES),
            "FC/without_dups_gene_level/counts.txt.summary"
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
    rule bbduk:
        input:
            r1 = "rawreads/{sample}_1.fq.gz",
            r2 = "rawreads/{sample}_2.fq.gz",
            fa = "FGS/rRNA.fa"
        output:
            r1 = "rawreads_rRNA_removed/{sample}_1.fq.gz",
            r2 = "rawreads_rRNA_removed/{sample}_2.fq.gz",
            stats = "rawreads_rRNA_removed/{sample}.stats.txt"
        params:
            threads = config["threads_bbduk"],
            kmer_length = config["kmer_length"]
        shell:
            "bbduk.sh k={params.kmer_length} threads={params.threads} in={input.r1} in2={input.r2} "
            "out={output.r1} out2={output.r2} ref={input.fa} stats={output.stats}"


    rule fastqc:
        input:
            expand("rawreads_rRNA_removed/{sample}_{paired}.fq.gz", sample=SAMPLES, paired=[1, 2])
        output:
            html="fastqc/raw/{sample}_{paired}_fastqc.html",
            zip="fastqc/raw/{sample}_{paired}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        params: "-t 2 --quiet"
        log:
            "logs/fastqc/raw/{sample}_{paired}.log"
        shell:
            "fastqc {params} {input} -o fastqc/raw/"


    rule trimmomatic:
        input:
            r1="rawreads_rRNA_removed/{sample}_1.fq.gz",
            r2="rawreads_rRNA_removed/{sample}_2.fq.gz"
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


    rule trimmed_fastqc:
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


    rule STAR:
        input:
            fq1 = "trimmed/{sample}.forward_paired.fq.gz",
            fq2 = "trimmed/{sample}.reverse_paired.fq.gz",
            dir = expand("FGS/{genome}", genome=config["genome"])
        output:
            # see STAR manual for additional output files -
            "star/{sample}.Aligned.sortedByCoord.out.bam"
        params:
            name = "star/{sample}.",
            threads = config["threads_star"],
            outReadsUnmapped = config["outReadsUnmapped"],
            outSAMattributes = config["outSAMattributes"],
            outSJfilterOverhangMin = config["outSJfilterOverhangMin"],
            outFilterMultimapNmax = config["outFilterMultimapNmax"],
            outFilterScoreMin = config["outFilterScoreMin"],
            outFilterMatchNminOverLread = config["outFilterMatchNminOverLread"],
            outFilterMismatchNmax = config["outFilterMismatchNmax"],
            outFilterMismatchNoverLmax = config["outFilterMismatchNoverLmax"],
            alignIntronMin = config["alignIntronMin"],
            alignIntronMax = config["alignIntronMax"],
            alignMatesGapMax = config["alignMatesGapMax"],
            alignSJoverhangMin = config["alignSJoverhangMin"],
            alignSJDBoverhangMin = config["alignSJDBoverhangMin"],
            alignSoftClipAtReferenceEnds = config["alignSoftClipAtReferenceEnds"],
            chimSegmentMin = config["chimSegmentMin"],
            chimScoreMin = config["chimScoreMin"],
            chimScoreSeparation = config["chimScoreSeparation"],
            chimJunctionOverhangMin = config["chimJunctionOverhangMin"],
            quantMode = config["quantMode"],
            twopassMode = config["twopassMode"],
            chimOutType = config["chimOutType"],
            outSAMtype = config["outSAMtype"],
            sjdbOverhang = config["sjdbOverhang"],
            outFilterMultimapScoreRange = config["outFilterMultimapScoreRange"]
            
            
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
            '--outReadsUnmapped {params.outReadsUnmapped}'
            '--outSAMattributes {params.outSAMattributes}'
            '--outSJfilterOverhangMin {params.outSJfilterOverhangMin}'
            '--outFilterMultimapNmax {params.outFilterMultimapNmax}'
            '--outFilterScoreMin {params.outFilterScoreMin}'
            '--outFilterMatchNminOverLread {params.outFilterMatchNminOverLread}'
            '--outFilterMismatchNmax {params.outFilterMismatchNmax}'
            '--outFilterMismatchNoverLmax {params.outFilterMismatchNoverLmax}'
            '--alignIntronMin {params.alignIntronMin}'
            '--alignIntronMax {params.alignIntronMax}'
            '--alignMatesGapMax {params.alignMatesGapMax}'
            '--alignSJoverhangMin {params.alignSJoverhangMin}'
            '--alignSJDBoverhangMin {params.alignSJDBoverhangMin}'
            '--alignSoftClipAtReferenceEnds {params.alignSoftClipAtReferenceEnds}'
            '--chimSegmentMin {params.chimSegmentMin}'
            '--chimScoreMin {params.chimScoreMin}'
            '--chimScoreSeparation {params.chimScoreSeparation}'
            '--chimJunctionOverhangMin {params.chimJunctionOverhangMin}'
            '--quantMode {params.quantMode}'
            '--twopassMode {params.twopassMode}'
            '--chimOutType {params.chimOutType}'
            '--outSAMtype = {params.outSAMtype}'
            '--sjdbOverhang = {params.sjdbOverhang}'
            '--outFilterMultimapScoreRange = {params.outFilterMultimapScoreRange}'


    rule multiqc:
        input:
            expand("fastqc/raw/{sample}_{paired}_fastqc.zip", sample=SAMPLES, paired=[1, 2]),
            expand("rawreads_rRNA_removed/{sample}_1.fq.gz", sample=SAMPLES),
            expand("rawreads_rRNA_removed/{sample}_2.fq.gz", sample=SAMPLES),
            expand("logs/trimmomatic/{sample}.log", sample=SAMPLES),
            expand("fastqc/trimmed/{sample}.{mate}_paired_fastqc.zip", sample=SAMPLES, mate=["forward", "reverse"]),
            expand("star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
            expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES),
            "FC/without_dups_gene_level/counts.txt.summary"
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
        "star/{sample}.Aligned.sortedByCoord.out.bam.csi"
    params:
        threads = config["threads_index_sorted_bams_with_dups"]
    shell:
        "samtools index -c -@ {params.threads} {input}"


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
        "removed_duplicates_alignments/{sample}.dedup.bam.csi"
    params:
        threads = config["threads_index_sorted_bams_without_dups"]
    shell:
        "samtools index -c -@ {params.threads} {input}"

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


rule featureCounts:
    input:
        with_dups_bams = lambda wildcards: expand("star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        without_dups_bams = expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES)
    output:
        with_dups_gene_level="FC/with_dups_gene_level/counts.txt",
        with_dups_exon_level="FC/with_dups_exon_level/counts.txt",
        without_dups_gene_level="FC/without_dups_gene_level/counts.txt",
        without_dups_exon_level="FC/without_dups_exon_level/counts.txt",
        without_dups_gene_level_summary="FC/without_dups_gene_level/counts.txt.summary"
    params:
        gtf = expand("FGS/{annotation}.gtf", annotation=config["annotation"])
    log:
        with_dups_gene_level="logs/FC/with_dups/with_dups_featurecount_gene.log",
        with_dups_exon_level="logs/FC/with_dups/with_dups_featurecount_exon.log",
        without_dups_gene_level="logs/FC/with_dups/with_dups_featurecount_gene.log",
        without_dups_exon_level="logs/FC/with_dups/with_dups_featurecount_exon.log"
#config[xxxx]
    shell:
        "featureCounts -T 4 -O -M -t exon -g gene_id -a {params.gtf} -o {output.with_dups_gene_level} {input.with_dups_bams} 2> {log.with_dups_gene_level} && "
        "featureCounts -T 4 -O -M -f -t exon -g gene_id -a {params.gtf} -o {output.with_dups_exon_level} {input.with_dups_bams} 2> {log.with_dups_exon_level} && "
        "featureCounts -T 4 -O -t exon -g gene_id -a {params.gtf} -o {output.without_dups_gene_level} {input.without_dups_bams} 2> {log.without_dups_gene_level} && "
        "featureCounts -T 4 -O -f -t exon -g gene_id -a {params.gtf} -o {output.without_dups_exon_level} {input.without_dups_bams} 2> {log.without_dups_exon_level}"
