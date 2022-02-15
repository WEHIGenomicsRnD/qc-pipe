if process_from_bcl:

    rule bcl2fastq:
        input:
            samplesheet=config["sample_sheet"],
        output:
            fastq=expand(
                "results/bcl_output/{sample}_{lane}_{readend}_001.fastq.gz",
                sample=samples,
                lane=lanes,
                readend=READENDS,
            ),
        log:
            "logs/bcl2fastq.log",
        envmodules:
            "bcl2fastq/2.20.0",
        threads: cluster["mergelanes"]["threads"]
        resources:
            mem_mb=cluster["mergelanes"]["mem_mb"],
            runtime=cluster["mergelanes"]["runtime"],
        shell:
            """
            bcl2fastq \
                --runfolder-dir {config[bcl_input]} \
                --output-dir results/bcl_output \
                --sample-sheet {input} \
                {config[params][bcl2fastq]}
            """


rule mergelanes:
    input:
        fastq=[
            "{bcl_dir}/{sample}_{lane}_{readend}_001.fastq.gz".format(
                bcl_dir=bcl_dir, sample=sample, lane=lane, readend=readend
            )
            for sample, lane, readend in zip(
                itertools.cycle(samples), itertools.cycle(lanes), READENDS
            )
        ],
    output:
        "fastq/{sample}_{readend}.fastq.gz",
    log:
        "logs/mergelanes_{sample}_{readend}.log",
    threads: cluster["mergelanes"]["threads"]
    resources:
        mem_mb=cluster["mergelanes"]["mem_mb"],
        runtime=cluster["mergelanes"]["runtime"],
    params:
        bcl_dir=bcl_dir,
    shell:
        "cat {params.bcl_dir}/{wildcards.sample}_*_{wildcards.readend}_001.fastq.gz > {output}"
