#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path
import sys


def load_bin_contig_map(tsv_file):
    """Load bin → contigs mapping from TSV file."""
    bin2contigs = defaultdict(set)

    with open(tsv_file) as f:
        header = next(f, None)  # skip header if present
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            bin_name, contig = parts[0], parts[1]
            bin2contigs[bin_name].add(contig)

    return bin2contigs


def split_fasta_by_bin(assembly_fasta, bin2contigs, outdir):
    """Stream assembly FASTA and write contigs into per-bin FASTA files."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    writers = {}
    contig_to_bins = defaultdict(list)

    # invert mapping for fast lookup
    for bin_name, contigs in bin2contigs.items():
        for c in contigs:
            contig_to_bins[c].append(bin_name)

    current_contig = None
    seq_buffer = []

    def flush():
        if current_contig is None:
            return
        if current_contig in contig_to_bins:
            for bin_name in contig_to_bins[current_contig]:
                if bin_name not in writers:
                    writers[bin_name] = open(
                        outdir / f"{bin_name}.fasta", "w"
                    )
                writers[bin_name].write(f">{current_contig}\n")
                writers[bin_name].write("".join(seq_buffer))

    with open(assembly_fasta) as f:
        for line in f:
            if line.startswith(">"):
                flush()
                current_contig = line[1:].strip().split()[0]
                seq_buffer = []
            else:
                seq_buffer.append(line)
        flush()

    for fh in writers.values():
        fh.close()


def main():
    parser = argparse.ArgumentParser(
        description="Split assembly FASTA into per-bin FASTA files"
    )
    parser.add_argument(
        "--assembly", required=True,
        help="Assembly FASTA file"
    )
    parser.add_argument(
        "--map", required=True,
        help="TSV file with columns: bin contig"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory for per-bin FASTA files"
    )

    args = parser.parse_args()

    if not Path(args.assembly).exists():
        sys.exit(f"ERROR: Assembly not found: {args.assembly}")
    if not Path(args.map).exists():
        sys.exit(f"ERROR: Mapping file not found: {args.map}")

    bin2contigs = load_bin_contig_map(args.map)
    if not bin2contigs:
        sys.exit("ERROR: No bin–contig mappings found")

    split_fasta_by_bin(
        assembly_fasta=args.assembly,
        bin2contigs=bin2contigs,
        outdir=args.outdir
    )


if __name__ == "__main__":
    main()
