if demultiplex:

    if demux_tool == "bcl2fastq":

        rule bcl2fastq:
            input:
                samplesheet=config["sample_sheet"],
            output:
                fastq=expand(
                    "results/demultiplexed/{sample}_{lane}_{readend}_001.fastq.gz",
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
                    --runfolder-dir {config[raw_input]} \
                    --output-dir results/demultiplexed \
                    --sample-sheet {input} \
                    {config[params][bcl2fastq]}
                """


    elif demux_tool == "mgi_splitbarcode":

        rule splitbarcodes:
            input:
                fq1="{raw_input}/{lane}/{sample_prefix}{lane}_read_1.fq.gz".format(
                    raw_input=config["raw_input"],
                    lane="{lane}",
                    sample_prefix="{sample_prefix}",
                ),
                fq2="{raw_input}/{lane}/{sample_prefix}{lane}_read_2.fq.gz".format(
                    raw_input=config["raw_input"],
                    lane="{lane}",
                    sample_prefix="{sample_prefix}",
                ),
                barcodes=config["params"]["mgi"]["bc_file"],
            output:
                fastq=expand(
                    "results/demultiplexed/{lane}/{sample_prefix}_{lane}_{sample_id}_{readend}.fq.gz",
                    lane="{lane}",
                    sample_prefix="{sample_prefix}",
                    sample_id=sample_ids,
                    readend=[1, 2],
                ),
            log:
                "logs/splitbarcode_{sample_prefix}_{lane}.log",
            threads: cluster["bcl2fastq"]["threads"]
            resources:
                mem_mb=cluster["bcl2fastq"]["mem_mb"],
                runtime=cluster["bcl2fastq"]["runtime"],
            params:
                binary=config["params"]["mgi"]["binary"],
                bs1=config["params"]["mgi"]["bs1"],
                bs2=config["params"]["mgi"]["bs2"],
                extra=config["params"]["mgi"]["extra"],
            shell:
                """
                {params.binary} \
                    -1 {input.fq1} \
                    -2 {input.fq2} \
                    -B {input.barcodes} \
                    -b {params.bs1} \
                    -b {params.bs2} \
                    -o results/demultiplexed/{wildcards.lane} \
                    -n {threads} \
                    {params.extra}
                """


    elif demux_tool == "bcl-convert":

        rule bcl2fastq:
            input:
                samplesheet=config["sample_sheet"],
            output:
                fastq=expand(
                    "results/demultiplexed/{sample}_{lane}_{readend}_001.fastq.gz",
                    sample=samples,
                    lane=lanes,
                    readend=READENDS,
                ),
            log:
                "logs/bclconvert.log",
            envmodules:
                "bcl-convert/3.9.3",
            threads: cluster["bcl2fastq"]["threads"]
            resources:
                mem_mb=cluster["bcl2fastq"]["mem_mb"],
                runtime=cluster["bcl2fastq"]["runtime"],
            shell:
                """
                bcl-convert \
                    --bcl-input-directory {config[raw_input]} \
                    --output-directory results/demultiplexed \
                    --sample-sheet {input} \
                    --force
                """


if demux_tool in ["bcl2fastq", "bcl-convert"]:

    rule mergelanes:
        input:
            fastqs=expand(
                "{demux_dir}/{sample}_{lane}_{readend}_001.fastq.gz",
                demux_dir=demux_dir,
                sample=samples,
                lane=lanes,
                readend=READENDS,
            ),
        output:
            "fastq/{sample}_{readend}.fastq.gz",
        log:
            "logs/mergelanes_{sample}_{readend}.log",
        threads: cluster["mergelanes"]["threads"]
        resources:
            mem_mb=cluster["mergelanes"]["mem_mb"],
            runtime=cluster["mergelanes"]["runtime"],
        params:
            demux_dir=demux_dir,
        shell:
            "cat {params.demux_dir}/{wildcards.sample}_*_{wildcards.readend}_001.fastq.gz > {output}"


elif demux_tool == "mgi_splitbarcode":

    rule mergelanes:
        input:
            fastqs=expand(
                "{demux_dir}/{lane}/{sample_prefix}_{lane}_{sample_id}_{readend}.fq.gz",
                demux_dir=demux_dir,
                lane=lanes,
                sample_prefix=sample_prefix,
                sample_id=sample_ids,
                readend=[1, 2],
            ),
        output:
            "fastq/{sample}_{readend}.fastq.gz",
        log:
            "logs/mergelanes_{sample}_{readend}.log",
        threads: cluster["mergelanes"]["threads"]
        resources:
            mem_mb=cluster["mergelanes"]["mem_mb"],
            runtime=cluster["mergelanes"]["runtime"],
        params:
            demux_dir=demux_dir,
            sample_prefix=sample_prefix,
            readend=lambda w: w.readend.split("R")[-1],
            sample_num=lambda w: w.sample.split("_")[-1],
        shell:
            """
            cat {params.demux_dir}/*/{params.sample_prefix}_*_{params.sample_num}_{params.readend}.fq.gz > {output}
            """
