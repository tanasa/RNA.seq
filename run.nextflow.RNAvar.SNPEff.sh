#!/bin/bash

# https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/

export CONDA_ENVS_PATH=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_envs
export CONDA_PKGS_DIRS=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_pkgs
export NXF_CONDA_CACHEDIR=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_cache
export TMPDIR=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_tmp
export TEMP=$TMPDIR
export TMP=$TMPDIR

nextflow run nf-core/rnavar \
-r 1.2.3 \
-profile apptainer \
--input RNAseq_Jaz_sample1_all_reads.rnavar.csv \
--outdir RNAseq_Jaz_sample1_all_reads.rnavar.results \
--genome GRCh38 \
--aligner "star" \
--star_index /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/STARIndex_2.7.11b \
--fasta /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa \
--dict /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.dict \
--fasta_fai /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa.fai \
--gtf /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genes.gtf \
-work-dir ./RNAseq_Jaz_sample1_all_reads.rnavar.work \
--remove_duplicates true \
-resume \
--dbsnp /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
--dbsnp_tbi /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
--known_indels /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known_indels_tbi /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
--snpeff_db GRCh38.99 \
--snpeff_cache /mnt/nfs/CX000008_DS1/projects/btanasa/snpeff_cache/ \
--tools "snpeff"

