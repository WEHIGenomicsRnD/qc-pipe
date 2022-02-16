import pandas as pd
import os
import yaml
import itertools
import re
from glob import iglob


# ------------- load cluster config ------------
with open("config/cluster.yaml", "r") as stream:
    try:
        cluster = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc, file=sys.stderr)

# ------------- globals ------------
READENDS = ["R1", "R2"]

# ------------- set up samples ------------
process_from_bcl = bool(config["process_from_bcl"])

bcl_dir = config["merge_from_dir"]
merge_without_demux = bcl_dir != "" and bcl_dir is not None
bcl_dir = bcl_dir if merge_without_demux else "results/bcl_output"

lanes = config["lanes"]
if lanes:
    lanes = ["L%s" % str(lane).zfill(3) for lane in range(1, int(lanes) + 1)]

if process_from_bcl:
    # process sample sheet for bcl2fastq
    from sample_sheet import SampleSheet

    sample_sheet = SampleSheet(config["sample_sheet"])
    samples = [
        "%s_S%d" % (s.sample_name, idx + 1)
        for idx, s in enumerate(sample_sheet.samples)
    ]

elif merge_without_demux:
    # in this case, demultiplexed fastqs exist, but are not merged
    fastq_dir = config["merge_from_dir"]
    fqs = iglob(f"{fastq_dir}/*_R1_001.fastq.gz")

    # extract basename of full file path
    base = [os.path.basename(i) for i in fqs]

    # extract sample name
    samples = []
    for f in base:
        samples.append(re.split("_L[0-9]{3}", f)[0])

else:
    # in this case, we expect the fastq files to already exist
    fastq_dir = config["fastq_dir"]
    fqs = iglob(f"{fastq_dir}/*_R1.fastq.gz")

    # extract basename of full file path
    base = [os.path.basename(i) for i in fqs]

    # extract sample name
    samples = []
    for f in base:
        samples.append(f.split("_R1")[0])

# ------------- output functions ------------


def get_mergelanes_output():
    mergelanes_output = expand(
        "fastq/{sample}_{readend}.fastq.gz", sample=samples, readend=READENDS
    )
    return mergelanes_output


def get_fastqc_output():
    fastqc_output = expand(
        "results/fastQC/{sample}_{readends}_fastqc.html",
        sample=samples,
        readends=READENDS,
    )
    return fastqc_output


def get_fastqscreen_output():
    fastqscreen_output = expand(
        "results/fastqScreen/{sample}_R2_screen.html", sample=samples
    )
    return fastqscreen_output


def get_align_output():
    align_output = expand("results/aligned/{sample}.bam", sample=samples)
    return align_output
