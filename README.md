# 🧬 protein_translation
Translating genes to protein sequences using gene annotation files (GTFs) and the ESM model.

## Set up environment

#### Option A: from 'environment.yml' file

```bash
conda env create -f environment.yaml
conda activate protein-env
```

#### Option B: from 'requirements.txt' file
```bash
conda create -n protein-env python=3.10 -y
conda activate protein-env
pip install -r requirements.txt
```

#### Option C: manual installation
```bash
conda create -n protein_env python=3.10 -y
conda activate protein-env

pip install pandas numpy biopython pyfaidx esm
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
```
💡 For CPU-only systems, you can omit the '--index-url' flag

## Project setup

#### Step 1
Download the hg38.p14 reference genome FASTA file from NCBI:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

#### Step 2
Generate process the fasta file by splitting the chromosomes (only want the main ones):
```bash
python 0_split_chroms.py /path/to/GCF_000001405.40_GRCh38.p14_genomic.fna
```

#### Step 3
Unzip the processed GTF files:
```bash
unzip -n 'resources/*.zip' -d resources/
```

## Scripts


change the file paths to resources 
cahnge the slurm output files so its not in the main script