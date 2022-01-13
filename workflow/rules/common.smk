import pandas as pd
import os
import yaml
from glob import iglob


#------------- load cluster config ------------
with open('config/cluster.yaml', 'r') as stream:
    try:
        cluster = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc, file=sys.stderr)

#------------- globals ------------
READENDS = ['R1', 'R2']

#------------- set up samples ------------
process_from_bcl = bool(config['process_from_bcl'])

if process_from_bcl:
    # process sample sheet for bcl2fastq
    from sample_sheet import SampleSheet

    sample_sheet = SampleSheet(config['sample_sheet'])
    samples = ['%s_S%d' % (s.sample_name, idx + 1) for idx, s in enumerate(sample_sheet.samples)]

    lanes = int(config['lanes'])
    lanes = ['L%s' % str(lane).zfill(3) for lane in range(1, lanes + 1)]

else:
    # in this case, we expect the fastq files to already exist
    fqs = iglob('fastq/*_R1.fastq.gz')

    # extract basename of full file path
    base = [os.path.basename(i) for i in fqs]

    # extract sample name
    samples = []
    for f in base:
        samples.append(f.split('_R1')[0])

#------------- output functions ------------
def get_bcl2fastq_output():
    bcl2fastq_output = expand(
            'results/bcl_output/{sample}_{lane}_{readend}_001.fastq.gz',
            sample=samples,
            lane=lanes,
            readend=READENDS
    )
    return bcl2fastq_output

def get_mergelanes_output():
    mergelanes_output = expand(
            'fastq/{sample}_{readend}.fastq.gz',
            sample=samples,
            readend=READENDS
        )
    return mergelanes_output

def get_fastqc_output():
    fastqc_output = expand(
        'results/fastQC/{sample}_{readends}_fastqc.html',
        sample=samples,
        readends=READENDS
    )
    return fastqc_output

def get_fastqscreen_output():
    fastqscreen_output = expand(
        'results/fastqScreen/{sample}_R2_screen.html',
         sample=samples
    )
    return fastqscreen_output
