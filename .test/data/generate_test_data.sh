#!/bin/sh

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/166/975/GCF_013166975.1_ASM1316697v1/GCF_013166975.1_ASM1316697v1_genomic.fna.gz \
#	&& gunzip GCF_013166975.1_ASM1316697v1_genomic.fna.gz
refname=GCF_013166975.1_ASM1316697v1_genomic

# generated using ART v2.5.8
# conda create -n art -c conda-forge -c bioconda art=2016.06.05
# conda activate art
art_illumina -ss HS25 -i ${refname}.fna -p -l 150 -f 2 -m 200 -s 10 -na -rs 1642032278 -o ecoli_R \
	&& gzip `ls *fq`

# change extension to fastq.gz
for fq in `ls *gz`; do sample=${fq%.fq.gz}; mv $fq ${sample}.fastq.gz ; done

mkdir -p ../fastq
mv *gz ../fastq

# build bowtie2 index for fastq screen
# conda deactivate
# conda create -n bowtie -c conda-forge -c bioconda bowtie2=2.4.4
# conda activate bowtie

bowtie2-build ${refname}.fna $refname 

# build bwa index for alignment
# conda deactivate
# conda create -n bwa -c conda-forge -c bioconda bwa=0.7.17
# conda activate bwa
bwa index ${refname}.fna

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/166/975/GCF_013166975.1_ASM1316697v1/GCF_013166975.1_ASM1316697v1_genomic.gff.gz && gunzip GCF_013166975.1_ASM1316697v1_genomic.gff.gz

# build star index for alignment
# conda deactivate
# conda create -n star -c conda-forge -c bioconda star=2.7.8a
# conda activate star
STAR --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles ${refname}.fna \
     --genomeSAindexNbases 10
