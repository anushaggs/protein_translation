#!/bin/bash
#SBATCH --job-name=extract_var_sequences         # create a short name for your job
#SBATCH --output=/mnt/home/aaggarwal/gates_proj/protein_translation/slurm_output/slurm-%A_%a.out  # stdout file
#SBATCH --nodes=1                 # request one node
#SBATCH --time=15:00         # set the maximum runtime (1 day)
#SBATCH --mem=512G                 # request 16GB of memory (adjust as needed)
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-user=aa8417@princeton.edu

python 2_extract_var_sequences.py \
  --gtf_transcript_path resources/hg38.p14.ncbiRefSeq.transcript_final_gtf.csv \
  --gtf_cds_path resources/hg38.p14.ncbiRefSeq.CDS_final_gtf.csv \
  --vcf_path /path/to/variants.vcf \
  --genome_fasta_dir /path/to/fasta_dir \
  --output_esm_fasta /path/to/output_variant_sequences.fasta
