# 3_extract_var_sequences.py
import torch
import pickle
from typing import List, Optional
import pandas as pd
from Bio.Seq import Seq
from pyfaidx import Fasta
import os

# Runs on CPU

def load_gtf(gtf_path):
    # load the GTF as pandas dataframe

    col_names = [
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes", "gene_id", "transcript_id", "entrez_id"
    ]
    
    gtf_df = pd.read_csv(gtf_path)

    return gtf_df

def load_vcf(vcf_path):
    # load the VCF file
    # works with ClinVar's VCF format

    with open(vcf_path, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith("#CHROM"):
                header_line = line.strip().lstrip("#").split("\t")
                break

    # read using the header line as column names
    vcf_df = pd.read_csv(
        vcf_path,
        sep="\t",
        comment="#",
        names=header_line,
        dtype={"CHROM": str, "POS": int, "REF": str, "ALT": str}
    )

    return vcf_df.head(50)

def extract_gene_regions(gtf_df):
    # filter to protein-coding genes and extract
    # - entrez id
    # - chrom
    # - tss
    # - gene start
    # - gene end

    entrez_id = gtf_df['entrez_id'].astype(str)
    chrom = gtf_df['seqname']
    strand = gtf_df['strand']
    tss = gtf_df['start'].where(strand == '+', gtf_df['end'])
    gene_start = gtf_df['start']
    gene_end = gtf_df['end']
    transcript_id = gtf_df['transcript_id']
    
    return zip(entrez_id, transcript_id, chrom, strand, tss, gene_start, gene_end)

def apply_variants_to_exon_seq(entrez_id, seq, seq_start, variants):
    # applies single bp substitutions (variants) to an exon sequnce

    seq = list(seq)
    for pos, ref, alt, in sorted(variants, key=lambda x: x[0]):
        rel_pos = pos - seq_start

        if rel_pos < 0 or rel_pos >= len(seq):
            print(f"[{entrez_id}] Variant {pos} not within exon window ({seq_start}-{seq_start + len(seq)})", flush=True)
            continue

        if len(ref) == 1 and len(alt) == 1:
            if seq[rel_pos].upper() != ref.upper():
                print(f"[{entrez_id}] Reference mismatch at position {pos}: expected {ref}, found {seq[rel_pos]}", flush=True)
                continue
            seq[rel_pos] = alt.upper()

    return ''.join(seq)

def extract_cds_sequence(transcript_id, chrom, strand, cds_df, fasta_index, variants=[]):
    transcript_df = cds_df[cds_df['transcript_id'] == transcript_id]
    cds_df = transcript_df[transcript_df['feature'] == 'CDS']
    start_df = transcript_df[transcript_df['feature'] == 'start_codon']
    stop_df = transcript_df[transcript_df['feature'] == 'stop_codon']

    if cds_df.empty or start_df.empty or stop_df.empty:
        print(f"[{transcript_id}] Missing CDS or codon annotations", flush=True)
        return None 

    if chrom not in fasta_index:
        print(f"[{transcript_id}] Chromosome {chrom} not found in FASTA", flush=True)
        return None

    genome = fasta_index[chrom]
    chrom_len = len(genome[chrom])

    start_codon_start = start_df.iloc[0]['start']
    start_codon_end = start_df.iloc[0]['end']
    stop_codon_start = stop_df.iloc[0]['start']
    stop_codon_end = stop_df.iloc[0]['end']

    if strand == '+':
        cds_df = cds_df[
            (cds_df['end'] >= start_codon_start) &
            (cds_df['start'] <= stop_codon_end)
        ].sort_values(by='start')
    else:
        cds_df = cds_df[
            (cds_df['end'] >= stop_codon_start) &
            (cds_df['start'] <= start_codon_end)
        ].sort_values(by='start', ascending=False)

    if cds_df.empty:
        print(f"[{transcript_id}] No CDSs after trimming by start/stop codons", flush=True)
        return None

    cds_rows = cds_df.to_dict('records') # to preserve the order
    seq_parts = []

    for i, row in enumerate(cds_rows):
        start = row['start']
        end = row['end']
        frame = int(row['frame']) if pd.notna(row['frame']) and str(row['frame']).isdigit() else 0

        if strand == '+':
            if start <= start_codon_start <= end:
                start = start_codon_start
            if start <= stop_codon_end <= end:
                end = stop_codon_end
        else:
            if start <= start_codon_end <= end:
                end = start_codon_end
            if start <= stop_codon_start <= end:
                start = stop_codon_start

        if end <= start:
            print(f"[{transcript_id}] Skipping — invalid CDS coordinates ({start}, {end})", flush=True)
            return None

        exon_seq = genome[chrom][start - 1:end].seq.upper()
        exon_seq = apply_variants_to_exon_seq(transcript_id, exon_seq, start, variants)

        if i == 0:
            exon_seq = exon_seq[frame:]
        
        #exon_seq = exon_seq[frame:]
        seq_parts.append(exon_seq)
    
    if strand == '-':
        seq_parts = seq_parts[::-1]
        full_seq = ''.join(seq_parts)
        full_seq = str(Seq(full_seq).reverse_complement())
    else:
        full_seq = ''.join(seq_parts)

    if not full_seq:
        print(f"[{transcript_id}] Sequence empty after CDS extraction", flush=True)
        return None

    return full_seq

def annotate_vcf_with_entrez(vcf_df, gtf_df):
    # for each variant in the vcf file, find which gene it falls in using the gtf file
    # returns a new VCF file with "entrez_id" column

    vcf_df = vcf_df.copy()
    vcf_df['entrez_id'] = None

    gene_regions = extract_gene_regions(gtf_df)

    for entrez_id, _, chrom, _, _, gene_start, gene_end in gene_regions:
        chrom = str(chrom)
        mask = (
            (vcf_df['Chrom'] == chrom) &
            (vcf_df['Pos'] >= gene_start) & 
            (vcf_df['Pos'] <= gene_end)
        )
        vcf_df.loc[mask, 'entrez_id'] = entrez_id

    return vcf_df[vcf_df['entrez_id'].notnull()].copy()

def build_variant_dict_by_entrez(annotated_vcf):
    # returns a dict: entrez_id --> list of (pos, ref, alt)
    # only applies variants for now, no indels (single bp substitution)

    variant_dict = {}
    for _, row in annotated_vcf.iterrows():
        eid = str(row['entrez_id'])
        pos = int(row['Pos'])
        ref = row['Ref']
        alt = row['Alt']
        if pd.isnull(ref) or pd.isnull(alt):
            continue

        # only keep variants (single bp substitutions) -- no indels
        if len(ref) == 1 and len(alt) == 1 and ref != '-' and alt != '-':
            variant_dict.setdefault(eid, []).append((pos, ref, alt))
            
    return variant_dict

def translate_to_protein(entrez_id, genic_sequence, warning=True):
    seq = genic_sequence.replace('\n', '').replace(' ', '').strip().upper()

    # we only want these warnings if the sequence is WT sequence, otherwise these warnings should be bypassed (variant sequence)
    if len(seq) % 3 != 0:
        if warning:
            print(f"[{entrez_id}] Sequence length not divisible by 3 — skipping.")
            return None

    protein = str(Seq(seq).translate(to_stop=False))

    if not protein.startswith('M'):
        if warning:
            print(f"[{entrez_id}] Protein does not start with 'M' — skipping.")
            return None

    if '*' in protein[:-1]:  # allow stop only at very end
        if warning:
            print(f"[{entrez_id}] Internal stop codon found — skipping.")
            return None

    return protein

def main(
    gtf_transcript_path: str,
    gtf_cds_path: str,
    vcf_path: str,
    genome_fasta_dir: str,
    output_esm_fasta: str,
):
    # load the genome
    fasta_files = {
        f.replace(".fa", ""): os.path.join(genome_fasta_dir, f)
        for f in os.listdir(genome_fasta_dir) if f.endswith(".fa")
    }
    fasta_index = {k: Fasta(v) for k, v in fasta_files.items()}

    # load GTF
    gtf_transcripts = load_gtf(gtf_transcript_path)
    gtf_cds = load_gtf(gtf_cds_path)

    # get gene regions for transcripts
    gene_regions = extract_gene_regions(gtf_transcripts)

    # load VCF
    vcf = load_vcf(vcf_path)

    # filer the vcf to variants/indels in gene regions and organize the variants per gene
    annotated_vcf = annotate_vcf_with_entrez(vcf, gtf_transcripts)
    variant_dict = build_variant_dict_by_entrez(annotated_vcf)

    # output file handles
    esm_out = open(output_esm_fasta, "w")

    # iterate through genes and get genomic & protein sequences
    for i, (entrez_id, transcript_id, chrom, strand, tss, gene_start, gene_end) in enumerate(gene_regions):
        
        # progress update
        if i % 1000 == 0 and i > 0:
            print(f"Processed {i} genes/proteins...", flush=True)
            
        # gene sequence for translation for ESM input
        wt_genic_sequence = extract_cds_sequence(transcript_id, chrom, strand, gtf_cds, fasta_index) # no variants

        if (wt_genic_sequence is None) and (entrez_id in variant_dict):
            print(f"Variant exists for gene {entrez_id}, but couldn't add it because WT sequence is incorrect", flush=True)
            continue
        elif (wt_genic_sequence is not None) and (entrez_id in variant_dict):
            variants = variant_dict[entrez_id]
            var_genic_sequence = extract_cds_sequence(transcript_id, chrom, strand, gtf_cds, fasta_index, variants=variants)
            var_protein_sequence = translate_to_protein(entrez_id, var_genic_sequence, warning=False)
            if var_protein_sequence:
                variant_strings = [f"{chrom}:g.{pos}{ref}>{alt}" for (pos, ref, alt) in variants]
                variant_label = ";".join(variant_strings)
                esm_out.write(f">{entrez_id}|{variant_label}\n{var_protein_sequence}\n")
        elif (wt_genic_sequence is not None) and (entrez_id not in variant_dict):
            wt_protein_sequence = translate_to_protein(entrez_id, wt_genic_sequence, warning=True)
            if wt_protein_sequence:
               esm_out.write(f">{entrez_id}\n{wt_protein_sequence}\n")
        else:
            continue

    esm_out.close()

    print("Fasta files written for sequence input into esm!")

if __name__ == "__main__":
    gtf_transcript_path = "/mnt/home/aaggarwal/ceph/gates_proj/ncbi_genome_hg38.p14/hg38.p14.ncbiRefSeq.transcript_final_gtf.csv"
    gtf_cds_path = "/mnt/home/aaggarwal/ceph/gates_proj/ncbi_genome_hg38.p14/hg38.p14.ncbiRefSeq.CDS_final_gtf.csv"
    genome_fasta_dir = "/mnt/home/aaggarwal/ceph/gates_proj/ncbi_genome_hg38.p14"
    output_esm_fasta = "/mnt/home/aaggarwal/ceph/gates_proj/GNN_model_files/fasta_files/esm_sequences.fasta"

    main(
        gtf_transcript_path,
        gtf_cds_path,
        genome_fasta_dir,
        output_esm_fasta,
    )