#!/bin/bash
#SBATCH --job-name=get_esm_embed    # create a short name for your job
#SBATCH --output=//mnt/home/aaggarwal/gates_proj/protein_translation/slurm_output/slurm-%A_%a.out
#SBATCH --gres=gpu:1  # request 4 GPUs
#SBATCH --constraint=a100-80gb&sxm4
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1  # increase number of CPUs per task
#SBATCH --time=2:00:00
#SBATCH --mem=512G
#SBATCH --mail-type=begin        # send mail when job begins
#SBATCH --mail-type=end          # send mail when job ends
#SBATCH --mail-type=fail         # send mail if job fails
#SBATCH --mail-user=aa8417@princeton.edu

module purge

conda activate target-genes

python 3_generate_esm_embeddings.py \
  --input_fasta output/esm_ref_sequences.fasta \
  --output_pkl output/esm_ref_embeddings.fasta 