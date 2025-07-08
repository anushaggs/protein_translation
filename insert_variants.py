import numpy as np
import os
import pandas as pd
from pyfaidx import Fasta
import pickle
import torch
import matplotlib.pyplot as plt

FULL_LENGTH = 40000
SEQUENCE_WINDOW = FULL_LENGTH // 2
EXTRA_FLANK = 500

# healthy kidney line: 'ACH-001310'
chosen_sample = "ACH-001398" # cancer kidney

vcf = pd.read_csv('/mnt/home/aaggarwal/ceph/gates_proj/depmap/OmicsSomaticMutations.csv')

vcf_sample = vcf[vcf['ModelID'] == chosen_sample]

gtf = pd.read_csv('/mnt/home/aaggarwal/ceph/gates_proj/genome_data_hg38/final_gtf_with_seqs.csv')
FASTA_DIR = "/mnt/home/aaggarwal/ceph/gates_proj/genome_data_hg38"
OUTPUT_FASTA = f"{FASTA_DIR}/hg38_18579_tss_sequences_40kbp_MUTATED_{chosen_sample}.fa"

# get all the chromosome FASTA files
fasta_files = {
    f.replace(".fa", ""): os.path.join(FASTA_DIR, f)
    for f in os.listdir(FASTA_DIR) if f.endswith(".fa")
}
fasta_index = {k: Fasta(v) for k, v in fasta_files.items()}

# parse the vcf file
def build_variant_dict(df):
    variant_dict = {}
    for _, row in df.iterrows():
        chrom = row["Chrom"]
        pos = int(row["Pos"])
        ref = row["Ref"]
        alt = row["Alt"]
        if pd.isnull(ref) or pd.isnull(alt): continue
        variant_dict.setdefault(chrom, []).append((pos, ref, alt))
    return variant_dict

vcf_variants = build_variant_dict(vcf_sample)

# variant application logic
def apply_variants_to_sequence(seq, seq_start, variants, tss, final_length=40000):
    seq = list(seq)
    delta = 0
    offset = 0
    tss_index = tss - seq_start

    for pos, ref, alt in sorted(variants, key=lambda x: x[0]):
        rel_pos = pos -1 - seq_start + offset
        if rel_pos < 0 or rel_pos + len(ref) > len(seq):
            print(f"[SKIP] out of bounds!", flush=True)
            continue
        current = ''.join(seq[rel_pos:rel_pos+len(ref)])
        if current.upper() != ref.upper():
            print(f"[SKIP] base mismatch at {chrom}:{pos} — expected '{ref}', found '{current}' in genome", flush=True)
            continue
        else:
            print(f"[APPLY] {chrom}:{pos}  {ref} → {alt}", flush=True)
        seq[rel_pos:rel_pos+len(ref)] = list(alt)
        delta += len(alt) - len(ref)

        if pos < tss:
            tss_index += delta

        offset += delta

    mutated = ''.join(seq)

    # recenter the sequence around the TSS
    # [WARNING] this portion assumes all indels would be covered by the 1000bp buffer (extra flanking regions).
    start = tss_index - SEQUENCE_WINDOW
    end = tss_index + SEQUENCE_WINDOW
    return mutated[start:end]

# write the mutated sequence to FASTA
with open(OUTPUT_FASTA, "w") as out_f:
    for i, row in gtf.iterrows():
        if i % 1000 == 0:
            print(f"Processing gene {i + 1} of {len(gtf)}...")

        chrom = row['seqname']
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        gene_id = row['query']

        if chrom not in fasta_index:
            print(f"Skipping {gene_id} — no FASTA for {chrom}")
            continue

        genome = fasta_index[chrom]
        chrom_len = len(genome[chrom])
        tss = start if strand == "+" else end

        seq_start = tss - SEQUENCE_WINDOW - EXTRA_FLANK
        seq_end = tss + SEQUENCE_WINDOW + EXTRA_FLANK

        if seq_start < 0:
            seq_end += abs(seq_start)
            seq_start = 0
        if seq_end > chrom_len:
            shift = seq_end - chrom_len
            seq_start -= shift
            seq_end = chrom_len
            if seq_start < 0:
                print(f"Skipping {gene_id} — cannot fit full window.")
                continue

        try:
            seq = genome[chrom][seq_start:seq_end].seq.upper()
            if len(seq) != (FULL_LENGTH + 2*EXTRA_FLANK):
                print(f"Skipping {gene_id} — final length {len(seq)}")
                continue
        except Exception as e:
            print(f"Error with {gene_id} at {chrom}:{seq_start}-{seq_end}: {e}")
            continue

        variants = [
            (pos, ref, alt)
            for (pos, ref, alt) in vcf_variants.get(chrom, [])
            if seq_start <= pos < seq_end
        ]

        if variants:
            final_seq = apply_variants_to_sequence(
                seq, seq_start, variants, tss, final_length=FULL_LENGTH
            )
            if len(final_seq) != FULL_LENGTH:
                print(f"Warning: {gene_id} mutated sequence length is {len(final_seq)}")
                continue
        else:
            seq_start = tss - SEQUENCE_WINDOW
            seq_end = tss + SEQUENCE_WINDOW
    
            if seq_start < 0:
                seq_end += abs(seq_start)
                seq_start = 0
            if seq_end > chrom_len:
                shift = seq_end - chrom_len
                seq_start -= shift
                seq_end = chrom_len
                if seq_start < 0:
                    print(f"Skipping {gene_id} — cannot fit full window.")
                    continue

            final_seq = genome[chrom][seq_start:seq_end].seq.upper()

            if len(final_seq) != FULL_LENGTH:
                    print(f"Warning: {gene_id} mutated sequence length is {len(final_seq)}")
                    continue
            
        out_f.write(f">{gene_id}\n{final_seq}\n")

print('all done', flush=True)

