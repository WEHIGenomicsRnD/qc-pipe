import pandas as pd
import os
import yaml
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
aligner = str(config["aligner"]).lower()
demultiplex = bool(config["demultiplex"])
demux_tool = config["demux_tool"].lower() if config["demux_tool"] else None

demux_dir = config["merge_from_dir"]
merge_without_demux = demux_dir != "" and demux_dir is not None
demux_dir = demux_dir if merge_without_demux else "results/demultiplexed"

lanes = config["lanes"]
if lanes:
    lanes = ["L%s" % str(lane).zfill(3) for lane in range(1, int(lanes) + 1)]

if demultiplex:

    if demux_tool == "bcl2fastq":
        # process sample sheet for bcl2fastq
        from sample_sheet import SampleSheet

        sample_sheet = SampleSheet(config["sample_sheet"])
        samples = [
            "%s_S%d" % (s.sample_name, idx + 1)
            for idx, s in enumerate(sample_sheet.samples)
        ]

    elif demux_tool == "bcl-convert":
        # we have to parse this by hand because bcl-convert uses
        # v2 sample sheets, which the python module does not support
        sample_info = []
        with open(config["sample_sheet"], "r") as ss:
            found_sample_info = False
            lines = ss.readlines()
            sample_info_idx = np.where(np.array(lines) == "[BCLConvert_Data]\n")[0]

            if len(sample_info_idx) == 0:
                print(
                    "Could not find [BCLConvert_Data] for demultiplexing. Check your sample sheet.",
                    file=sys.stderr,
                )
                sys.exit()

            for line in lines[sample_info_idx[0] + 1 :]:
                if line == "\n":
                    break
                sample_info.append(line.rstrip())

            sample_info = [s.split(",") for s in sample_info]
            sample_info = pd.DataFrame(sample_info[1:], columns=sample_info[0])
            samples = [
                "%s_S%d" % (s, idx + 1) for idx, s in enumerate(sample_info.Sample_ID)
            ]

    elif demux_tool == "mgi_splitbarcode":

        raw_input = config["raw_input"]
        fqs = iglob(f"{raw_input}/*/*_1.fq.gz")

        # extract basename of full file path
        sample_prefix = [os.path.basename(i) for i in fqs]
        sample_prefix = re.split("_L[0-9]{2}", sample_prefix[0])[0]

        # infer lanes
        lanes = iglob(f"{raw_input}/*")
        lanes = [os.path.basename(lane) for lane in lanes]

        # make sure this includes only lane names
        assert all([bool(re.match("L[0-9]{2}", lane)) for lane in lanes])

        # construct sample names from base name + barcodes file
        barcodes = pd.read_csv(
            config["params"]["mgi"]["bc_file"], sep="\t", header=None, dtype="str"
        )
        sample_ids = barcodes[0].values
        samples = [f"{sample_prefix}_{sample_id}" for sample_id in sample_ids]

    else:
        print(
            f"ERROR: invalid demultiplexing tool specified: {demux_tool}",
            file=sys.stderr,
        )

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


# ------------- input/output functions ------------


def get_splitbarcodes_output():
    splitbarcodes_output = (
        expand(
            "results/demultiplexed/{lane}/{sample_prefix}_{lane}_{sample_id}_{readend}.fq.gz",
            lane=lanes,
            sample_prefix=sample_prefix,
            sample_id=sample_ids,
            readend=[1, 2],
        ),
    )
    return splitbarcodes_output


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
