rule align:
    input:
        r1='fastq/{sample}_R1.fastq.gz',
        r2='fastq/{sample}_R2.fastq.gz',
    output:
        'results/aligned/{sample}.bam'
    log:
        'logs/align_{sample}.log'
    conda:
        '../envs/bwa_samtools.yaml'
    threads:
        cluster['align']['threads']
    resources:
        mem_mb=cluster['align']['mem_mb'],
        runtime=cluster['align']['runtime'],
    shell:
        '''
        bwa mem -p -t {threads} {config[ref]} {input} | samtools view -bS - > {output}
        '''
