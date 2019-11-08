def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    
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
