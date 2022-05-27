rule fastQC:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
    output:
        r1="results/fastQC/{sample}_R1_fastqc.html",
        r2="results/fastQC/{sample}_R2_fastqc.html",
    log:
        "logs/fastQC_{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    threads: cluster["fastqc"]["threads"]
    resources:
        mem_mb=cluster["fastqc"]["mem_mb"],
        runtime=cluster["fastqc"]["runtime"],
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -o results/fastQC
        """


rule fastqScreen:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
        config_file=config["FastqScreen_config"],
    output:
        r1="results/fastqScreen/{sample}_R1_screen.html",
        r2="results/fastqScreen/{sample}_R2_screen.html",
    log:
        "logs/fastqScreen_{sample}.log",
    conda:
        "../envs/fastqscreen.yaml"
    threads: cluster["fastqscreen"]["threads"]
    resources:
        mem_mb=cluster["fastqscreen"]["mem_mb"],
        runtime=cluster["fastqscreen"]["runtime"],
    shell:
        """
        fastq_screen --conf {input.config_file} \
            --outdir results/fastqScreen \
            {input.r1} {input.r2}
        """


if aligner != "none":

    rule samtools_stats:
        input:
            "results/aligned/{sample}.bam",
        output:
            stats="results/samtools/stats/{sample}.txt",
            flagstat="results/samtools/flagstat/{sample}.txt",
        log:
            "logs/samtools_stats_{sample}.log",
        conda:
            "../envs/bwa_samtools.yaml"
        threads: cluster["samtools_stats"]["threads"]
        resources:
            mem_mb=cluster["samtools_stats"]["mem_mb"],
            runtime=cluster["samtools_stats"]["runtime"],
        shell:
            """
            samtools stats {input} > {output.stats}
            samtools flagstat {input} > {output.flagstat}
            """

    rule qualimap:
        input:
            "results/aligned/{sample}.bam",
        output:
            "results/qualimap/{sample}/qualimapReport.html",
        log:
            "logs/qualimap_{sample}.log",
        conda:
            "../envs/qualimap.yaml"
        threads: cluster["qualimap"]["threads"]
        resources:
            mem_mb=cluster["qualimap"]["mem_mb"],
            runtime=cluster["qualimap"]["runtime"],
        shell:
            """
            qualimap bamqc \
                --java-mem-size={resources.mem_mb}M \
                -bam {input} \
                -outdir results/qualimap/{wildcards.sample} \
                -outformat HTML \
                -nt {threads}
            """

    rule multiQC:
        input:
            fastqc=expand("results/fastQC/{sample}_R1_fastqc.html", sample=samples),
            fastqscreen=expand(
                "results/fastqScreen/{sample}_R1_screen.html", sample=samples
            ),
            samtools=expand("results/samtools/stats/{sample}.txt", sample=samples),
            qualimap=expand(
                "results/qualimap/{sample}/qualimapReport.html", sample=samples
            ),
        output:
            "results/multiqc/multiqc_report.html",
        log:
            "logs/multiQC.log",
        conda:
            "../envs/multiqc.yaml"
        threads: cluster["multiqc"]["threads"]
        resources:
            mem_mb=cluster["multiqc"]["mem_mb"],
            runtime=cluster["multiqc"]["runtime"],
        shell:
            """
            dirs="results/fastQC/ \
                  results/fastqScreen/ \
                  results/qualimap/ \
                  results/samtools/"

            if [ -d results/bcl_output ]; then
                dirs="results/bcl_output $dirs"
            fi

            multiqc $dirs -o results/multiqc -f
            """


else:

    rule multiQC:
        input:
            fastqc=expand("results/fastQC/{sample}_R1_fastqc.html", sample=samples),
            fastqscreen=expand(
                "results/fastqScreen/{sample}_R1_screen.html", sample=samples
            ),
        output:
            "results/multiqc/multiqc_report.html",
        log:
            "logs/multiQC.log",
        conda:
            "../envs/multiqc.yaml"
        threads: cluster["multiqc"]["threads"]
        resources:
            mem_mb=cluster["multiqc"]["mem_mb"],
            runtime=cluster["multiqc"]["runtime"],
        shell:
            """
            dirs="results/fastQC/ \
                  results/fastqScreen/"

            if [ -d results/bcl_output ]; then
                dirs="results/bcl_output $dirs"
            fi

            multiqc $dirs -o results/multiqc -f
            """
