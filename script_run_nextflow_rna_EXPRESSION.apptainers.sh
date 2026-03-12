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

nextflow run nf-core/rnaseq \
-r 3.14.0 \
-profile apptainer \
--input RNAseq_Jaz.csv \
--outdir RNAseq_Jaz.results \
--genome GRCh38 \
--aligner "star_rsem" \
--pseudo_aligner "salmon" \
-work-dir ./RNAseq_Jaz.work \
-resume \
--remove_ribo_rna false \
--stringtie_ignore_gtf true \
--save_unaligned true \
--gtf  /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genes.gtf \
--star_index /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/STARIndex/ \
--fasta /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa \
--dict /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.dict \
--fasta_fai /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genome.fa.fai \
--gene_bed /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/genes.bed \
--rsem_index /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/RSEMIndex \
--salmon_index /mnt/nfs/CX000008_DS1/projects/btanasa/iGenomes/GRCh38/SalmonIndex \
--contaminant_screening true \
--skip_dupradar

# --skip_deseq_qc
# --skip_biotype_qc
# --skip_sortmerna

# --skip_rseqc 
# --contaminant_screening true
# -c /mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/nextflow/nf.config