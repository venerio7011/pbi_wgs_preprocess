#!/usr/bin/env python3
"""
GenBank Conversion Utilities

Convert GenBank files to various FASTA formats.
"""

import sys
from pathlib import Path
from typing import Optional
from Bio import SeqIO


def genbank_to_faa(genbank_file: str, faa_file: str) -> int:
    """
    Extract protein sequences from GenBank file and write to .faa format.

    Args:
        genbank_file: Path to input GenBank file
        faa_file: Path to output .faa file

    Returns:
        Number of protein sequences extracted
    """
    protein_count = 0

    with open(faa_file, 'w') as out_handle:
        for record in SeqIO.parse(genbank_file, 'genbank'):
            for feature in record.features:
                if feature.type == 'CDS':
                    if 'translation' in feature.qualifiers:
                        protein_seq = feature.qualifiers['translation'][0]

                        if 'locus_tag' in feature.qualifiers:
                            protein_id = feature.qualifiers['locus_tag'][0]
                        elif 'gene' in feature.qualifiers:
                            protein_id = feature.qualifiers['gene'][0]
                        elif 'protein_id' in feature.qualifiers:
                            protein_id = feature.qualifiers['protein_id'][0]
                        else:
                            protein_id = f"protein_{protein_count + 1}"

                        product = feature.qualifiers.get('product', ['hypothetical protein'])[0]

                        out_handle.write(f">{protein_id} {product}\n")
                        out_handle.write(f"{protein_seq}\n")
                        protein_count += 1

    return protein_count


def genbank_to_cds_fasta(genbank_file: str, output_file: str) -> int:
    """
    Extract CDS nucleotide sequences from GenBank file.

    Args:
        genbank_file: Path to input GenBank file
        output_file: Path to output FASTA file

    Returns:
        Number of CDS sequences extracted
    """
    cds_count = 0

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(genbank_file, 'genbank'):
            for feature in record.features:
                if feature.type == 'CDS':
                    if 'locus_tag' in feature.qualifiers:
                        cds_id = feature.qualifiers['locus_tag'][0]
                    elif 'gene' in feature.qualifiers:
                        cds_id = feature.qualifiers['gene'][0]
                    else:
                        cds_id = f"CDS_{cds_count + 1}"

                    product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
                    cds_seq = feature.extract(record.seq)

                    out_handle.write(f">{cds_id} {product}\n")
                    out_handle.write(f"{str(cds_seq)}\n")
                    cds_count += 1

    return cds_count


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python genbank_utils.py <input.genbank> <output.faa>")
        print("       python genbank_utils.py <input.genbank> <output.fna> --cds")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if len(sys.argv) > 3 and sys.argv[3] == '--cds':
        count = genbank_to_cds_fasta(input_file, output_file)
        print(f"Extracted {count} CDS nucleotide sequences to {output_file}")
    else:
        count = genbank_to_faa(input_file, output_file)
        print(f"Extracted {count} protein sequences to {output_file}")
