#!/bin/bash

# https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/

# --- Conda / Nextflow cache ---
# export CONDA_ENVS_PATH=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_envs
# export CONDA_PKGS_DIRS=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_pkgs
# export NXF_CONDA_CACHEDIR=/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/conda_cache

# --- Put temporary files on local scratch ---
export TMPDIR=/mnt/local/scratch/btanasa/nextflow_tmp
export TEMP=$TMPDIR
export TMP=$TMPDIR

# Optional but often helpful
export NXF_TEMP=$TMPDIR

# Create the temp directory if needed
mkdir -p "$TMPDIR"

# --- Apptainer / Singularity binds ---
export APPTAINER_BIND="/mnt/nfs/CX000008_DS1/projects/jaz/jharris3/data:/mnt/nfs/CX000008_DS1/projects/jaz/jharris3/data,/mnt/local/scratch:/mnt/local/scratch"
export SINGULARITY_BIND="$APPTAINER_BIND"
export NXF_APPTAINER_CACHEDIR=/mnt/nfs/CX000008_DS1/projects/btanasa/containers

nextflow run nf-core/rnavar \
-r 1.2.3 \
-profile apptainer \
--input RNAseq_Jaz.rnavar.csv \
--outdir RNAseq_Jaz.rnavar.results \
--genome GRCh38 \
--aligner "star" \
--star_index /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/STARIndex_2.7.11b \
--fasta /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa \
--dict /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.dict \
--fasta_fai /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa.fai \
--gtf /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genes.gtf \
-work-dir ./RNAseq_Jaz.rnavar.work \
--remove_duplicates true \
-resume \
--dbsnp /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
--dbsnp_tbi /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
--known_indels /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known_indels_tbi /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
--snpeff_db GRCh38.99 \
--snpeff_cache /mnt/nfs/CX000008_DS1/projects/btanasa/snpeff_cache/ \
--tools "snpeff"

# --snpeff_cache /mnt/nfs/CX000008_DS1/projects/btanasa/snpeff_cache/ \
# --vep_cache /mnt/nfs/CX000008_DS1/projects/btanasa/vep_cache/ 

# --snpeff_db GRCh38.105 \
# --vep_genome 102_GRCh38 \

# --vep_cache /mnt/nfs/CX000008_DS1/projects/btanasa/vep_cache \
# --snpeff_cache /mnt/nfs/CX000008_DS1/projects/btanasa/snpeff_cache/

# --vep_cache /mnt/nfs/CX000008_DS1/projects/btanasa/vep_cache/102_GRCh38/homo_sapiens/102_GRCh38 \
# --vep_dbnsfp true \
# --dbnsfp
# --dbnsfp_tbi
# --snpeff_cache /mnt/nfs/CX000008_DS1/projects/btanasa/snpeff_cache/GRCh38.105/GRCh38.105

# --extract_umi false
# --vep_loftee
# --spliceai_snv
# --spliceai_snv_tbi
# --spliceai_indel
# --spliceai_indel_tbi
# -c /mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/nf.config