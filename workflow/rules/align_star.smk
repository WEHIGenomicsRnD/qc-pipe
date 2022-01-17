rule align:
    input:
        fq1="fastq/{sample}_R1.fastq.gz",
        fq2="fastq/{sample}_R2.fastq.gz",
        gtf=config["gtf"],
    output:
        "results/aligned/{sample}.bam",
    log:
        "logs/align_{sample}.log",
    conda:
        "../envs/star.yaml"
    threads: cluster["align"]["threads"]
    resources:
        mem_mb=cluster["align"]["mem_mb"],
        runtime=cluster["align"]["runtime"],
    params:
        index=config["star_index"],
        extra=config["params"]["star"],
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand zcat \
            --outFileNamePrefix results/aligned/{wildcards.sample}. \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbGTFfile {input.gtf} \
            {params.extra}

        cd results/aligned && 
            ln -s {wildcards.sample}.Aligned.sortedByCoord.out.bam \
                {wildcards.sample}.bam
        """
