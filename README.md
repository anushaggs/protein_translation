# 🧬 protein_translation
Translating genes to protein sequences using gene annotation files (GTFs) and the ESM model.

## Set up environment (with GPU)

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

## Scripts overview

### `1_extract_ref_sequences.py`

**Description**  
Extracts reference coding sequences (CDS) for canonical transcripts and translates them into protein sequences using the hg38.p14 reference genome. CPU-only task.

**Inputs**
- Transcript-level GTF file (`*_transcript_final_gtf.csv`)
- CDS-level GTF file (`*_CDS_final_gtf.csv`)
- Genome FASTA directory (with per-chromosome `.fa` files)

**Output**
- FASTA file of reference protein sequences formatted for ESM input

**Usage**
```bash
python 1_extract_ref_sequences.py \
  --gtf_transcript_path resources/hg38.p14.ncbiRefSeq.transcript_final_gtf.csv \
  --gtf_cds_path resources/hg38.p14.ncbiRefSeq.CDS_final_gtf.csv \
  --genome_fasta_dir path/to/genome_fasta/ \
  --output_esm_fasta output/esm_ref_sequences.fasta
```

### `2_extract_var_sequences.py`

**Description**  
Extracts and translates coding sequences (CDS) from canonical transcripts in a GTF file using a genome FASTA reference. For transcripts that have variants (from a VCF file), it applies the variants and outputs translated variant protein sequences. For transcripts without any associated variants, it outputs their reference (wild-type) protein sequences. All outputs are formatted for ESM input. CPU-only task.

Variants are limited to single bp substitutions. Indels are not handeled.

**Inputs**
- Transcript-level GTF file (`*_transcript_final_gtf.csv`)
- CDS-level GTF file (`*_CDS_final_gtf.csv`)
- Genome FASTA directory (with per-chromosome `.fa` files)
- VCF file (`.vcf`file)

**Output**
- FASTA file of translated protein sequences:
  - Variant protein sequences for transcripts with mapped single-nucleotide variants (SNVs)
    - Variants are listed in the header (e.g., `>12345|chr7:g.117199644A>G`)
  - Reference protein sequences for transcripts with **no variants**

**Usage**
```bash
python 2_extract_var_sequences.py \
  --gtf_transcript_path resources/hg38.p14.ncbiRefSeq.transcript_final_gtf.csv \
  --gtf_cds_path resources/hg38.p14.ncbiRefSeq.CDS_final_gtf.csv \
  --vcf_path /path/to/variants.vcf \
  --genome_fasta_dir /path/to/fasta_dir \
  --output_esm_fasta /path/to/output_variant_sequences.fasta
```

### `3_generate_esm_embeddings.py`

**Description**  
Takes the reference or variant FASTA sequences and computes ESM model embeddings for downstream tasks. CPU or GPU task (GPU will automatically get detected by the script).

**Inputs**
- FASTA file from `1_extract_ref_sequences.py` or `2_extract_var_sequences.py`

**Output**
- Pickled dictionary of ESM embeddings (CLS token only)
- For protein sequences longer than 2048, the proteins are chunked and averaged

**Usage**
```bash
python 3_generate_esm_embeddings.py \
  --input_fasta output/esm_ref_sequences.fasta \
  --output_pkl output/esm_ref_embeddings.fasta 
```
