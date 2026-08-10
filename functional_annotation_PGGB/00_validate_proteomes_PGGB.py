#!/usr/bin/env python3
"""Validate and normalize one protein FASTA file."""

import argparse
import csv
import re
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("input_fasta")
parser.add_argument("output_fasta")
parser.add_argument("id_map")
parser.add_argument("statistics")
args = parser.parse_args()

ALLOWED_RESIDUES = set("ABCDEFGHIKLMNPQRSTVWXYZ*UOJBZ")


def read_fasta(path):
    header = None
    sequence = []
    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise ValueError("Sequence found before the first FASTA header")
                sequence.append(re.sub(r"\s+", "", line).upper())
    if header is not None:
        yield header, "".join(sequence)


try:
    Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.id_map).parent.mkdir(parents=True, exist_ok=True)
    Path(args.statistics).parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    lengths = []
    invalid_residues = 0

    with open(args.output_fasta, "w", newline="\n") as output, open(
        args.id_map, "w", newline="", encoding="utf-8"
    ) as mapping:
        writer = csv.writer(mapping, delimiter="\t")
        writer.writerow(["protein_id", "original_header"])

        for original_header, sequence in read_fasta(args.input_fasta):
            protein_id = original_header.split()[0]
            if not protein_id:
                raise ValueError("Empty FASTA identifier")
            if re.search(r"[^A-Za-z0-9_.:|+-]", protein_id):
                raise ValueError(f"Unsupported character in identifier: {protein_id}")
            if protein_id in seen:
                raise ValueError(f"Duplicated protein identifier: {protein_id}")
            if not sequence:
                raise ValueError(f"Empty sequence: {protein_id}")

            seen.add(protein_id)
            lengths.append(len(sequence))
            invalid_residues += sum(residue not in ALLOWED_RESIDUES for residue in sequence)
            output.write(f">{protein_id}\n")
            output.write("\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)) + "\n")
            writer.writerow([protein_id, original_header])

    if not lengths:
        raise ValueError("The input FASTA contains no sequences")
    if invalid_residues:
        raise ValueError(f"Found {invalid_residues} unsupported residues")

    total = sum(lengths)
    with open(args.statistics, "w", newline="", encoding="utf-8") as statistics:
        writer = csv.writer(statistics, delimiter="\t")
        writer.writerow(["proteins", "total_aa", "min_length", "mean_length", "max_length", "invalid_residues"])
        writer.writerow([len(lengths), total, min(lengths), f"{total / len(lengths):.2f}", max(lengths), invalid_residues])
except Exception as error:
    print(f"ERROR: {error}", file=sys.stderr)
    sys.exit(1)

