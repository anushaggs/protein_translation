import argparse
from pathlib import Path
from Bio import SeqIO

def load_mapping(mapping_path):
    mapping = {}
    with open(mapping_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                nc_id, chr_name = parts
                mapping[nc_id] = chr_name
    return mapping

def extract_chromosomes(fasta_path, mapping):
    output_dir = fasta_path.parent
    records_written = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        nc_id = record.id
        if nc_id in mapping:
            chr_name = mapping[nc_id]
            record.id = chr_name
            record.description = ""
            out_path = output_dir / f"{chr_name}.fa"
            with open(out_path, "w") as f:
                SeqIO.write(record, f, "fasta")
            records_written += 1

    print(f"Extracted {records_written} chromosomes to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Split FASTA into chromosome files using chrom_map")
    parser.add_argument(
        "fasta",
        help="Path to input FASTA file (e.g., GCF_000001405.40_GRCh38.p14_genomic.fna)"
    )

    args = parser.parse_args()
    fasta_path = Path(args.fasta)
    mapping_path = Path("resources/chrom_map.txt")

    if not mapping_path.exists():
        raise FileNotFoundError(f"Mapping file not found: {mapping_path}")

    mapping = load_mapping(mapping_path)
    extract_chromosomes(fasta_path, mapping)

if __name__ == "__main__":
    main()
