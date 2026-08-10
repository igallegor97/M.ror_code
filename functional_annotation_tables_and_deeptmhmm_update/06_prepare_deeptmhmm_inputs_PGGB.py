#!/usr/bin/env python3
"""Validate and split the five proteomes into traceable DeepTMHMM batches."""

import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--input-dir", required=True)
parser.add_argument("--output-dir", required=True)
parser.add_argument("--batch-size", type=int, default=2000)
args = parser.parse_args()

input_dir = Path(args.input_dir)
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
allowed = set("ACDEFGHIKLMNPQRSTVWY")
manifest_rows = []


def read_fasta(path):
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.upper())
    if header is not None:
        yield header, "".join(sequence)


for fasta in sorted(input_dir.glob("*.all_proteins.faa")):
    sample = fasta.name.removesuffix(".all_proteins.faa")
    records = list(read_fasta(fasta))
    invalid = [
        (identifier, sorted(set(sequence) - allowed))
        for identifier, sequence in records
        if set(sequence) - allowed
    ]
    if invalid:
        example = "; ".join(f"{identifier}:{''.join(chars)}" for identifier, chars in invalid[:10])
        raise SystemExit(
            f"ERROR: {sample} contains {len(invalid)} sequences with characters "
            f"not accepted by the official DeepTMHMM server. Examples: {example}"
        )

    sample_dir = output_dir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    for start in range(0, len(records), args.batch_size):
        batch_records = records[start : start + args.batch_size]
        batch_number = start // args.batch_size + 1
        batch_name = f"{sample}.deeptmhmm.batch_{batch_number:03d}.faa"
        batch_path = sample_dir / batch_name
        with batch_path.open("w", newline="\n", encoding="utf-8") as handle:
            for identifier, sequence in batch_records:
                handle.write(f">{identifier}\n")
                handle.write("\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)) + "\n")
        manifest_rows.append({
            "sample_id": sample,
            "batch_id": f"{batch_number:03d}",
            "input_fasta": str(batch_path),
            "proteins": len(batch_records),
            "amino_acids": sum(len(sequence) for _, sequence in batch_records),
            "expected_output": f"{sample}.batch_{batch_number:03d}.predicted_topologies.3line",
        })

manifest = output_dir / "deeptmhmm_batch_manifest.tsv"
with manifest.open("w", newline="", encoding="utf-8") as handle:
    fields = ["sample_id", "batch_id", "input_fasta", "proteins", "amino_acids", "expected_output"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

print(f"DeepTMHMM batches: {len(manifest_rows)}")
print(f"Proteins represented: {sum(row['proteins'] for row in manifest_rows)}")
print(f"Manifest: {manifest}")

