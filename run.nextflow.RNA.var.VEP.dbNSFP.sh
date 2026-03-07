#!/bin/bash

# https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/

export CONDA_ENVS_PATH=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_envs
export CONDA_PKGS_DIRS=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_pkgs
export TMPDIR=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_tmp
export TEMP=$TMPDIR
export TMP=$TMPDIR

nextflow run nf-core/rnavar \
-r 1.2.3 \
-profile apptainer \
--input RNAseq_Jaz_sample_500K.rnavar.csv \
--outdir RNAseq_Jaz_sample_500K.rnavar.results \
--genome GRCh38 \
--aligner "star" \
--star_index /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/STARIndex_2.7.11b \
--fasta /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa \
--dict /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.dict \
--fasta_fai /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa.fai \
--gtf /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genes.gtf \
-work-dir ./RNAseq_Jaz_sample_500K.rnavar.work  \
--remove_duplicates true \
-resume \
--dbsnp /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
--dbsnp_tbi /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
--known_indels /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known_indels_tbi /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
--tools "vep" \
--vep_genome GRCh38 \
--vep_cache_version 115 \
--vep_species homo_sapiens \
--vep_cache /mnt/nfs/CX000008_DS1/projects/btanasa/vep_cache/ \
--vep_custom_args "--everything --offline --format vcf" \
--vep_dbnsfp true \
--dbnsfp /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/dbNSFP/dbNSFP5.3.1a_grch38.gz \
--dbnsfp_tbi  /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/dbNSFP/dbNSFP5.3.1a_grch38.gz.tbi \
--dbnsfp_fields "rs_dbSNP,HGVSc_VEP,HGVSp_VEP,1000Gp3_EAS_AF,1000Gp3_AMR_AF,LRT_score,GERP++_RS,gnomAD_exomes_AF" 
