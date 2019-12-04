rule star_index:
    input:
        fasta = "FGS/{genome}.fasta"
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
        "0.40.2/bio/star/index"
