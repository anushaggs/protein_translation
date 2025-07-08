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

Step 1: Load the conda env

Step 2: Get the genome fasta file and separate into chromosomes

Step 3: unzip the gtf files in resources

Run 1, 2, 3, 

change the file paths to resources 
cahnge the slurm output files so its not in the main script