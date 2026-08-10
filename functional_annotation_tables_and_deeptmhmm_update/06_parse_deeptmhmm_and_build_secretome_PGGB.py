#!/usr/bin/env python3
"""Parse DeepTMHMM 3-line outputs and build a SignalP-filtered secretome."""

import argparse
import csv
import re
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--results", required=True)
parser.add_argument("--signalp-dir", required=True)
parser.add_argument("--proteome-dir", required=True)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()

deep_root = Path(args.results)
signalp_dir = Path(args.signalp_dir)
proteome_dir = Path(args.proteome_dir)
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    identifier = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence)
                identifier = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if identifier is not None:
        yield identifier, "".join(sequence)


def topology_segments(topology, state):
    segments = []
    for match in re.finditer(f"{re.escape(state)}+", topology):
        segments.append((match.start() + 1, match.end()))
    return segments


def parse_three_line(path):
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) % 3:
        raise ValueError(f"{path}: expected a multiple of three non-empty lines")
    for index in range(0, len(lines), 3):
        header, sequence, topology = lines[index : index + 3]
        if not header.startswith(">"):
            raise ValueError(f"{path}: invalid header at record {index // 3 + 1}")
        if len(sequence) != len(topology):
            raise ValueError(f"{path}: sequence/topology length mismatch for {header}")
        header_text = header[1:]
        protein_id = header_text.split()[0]
        prediction = header_text.split("|")[-1].strip() if "|" in header_text else ""
        yield protein_id, sequence, topology, prediction


deeptmhmm = {}
for path in sorted(deep_root.rglob("*.3line")):
    for protein_id, sequence, topology, prediction in parse_three_line(path):
        if protein_id in deeptmhmm:
            raise SystemExit(f"ERROR: duplicated DeepTMHMM protein ID: {protein_id}")
        tm_segments = topology_segments(topology, "M")
        deeptmhmm[protein_id] = {
            "sequence": sequence,
            "topology": topology,
            "prediction": prediction,
            "tm_segments": tm_segments,
            "source_file": str(path),
        }

if not deeptmhmm:
    raise SystemExit(f"ERROR: no *.3line files found under {deep_root}")

summary_fields = [
    "sample_id", "protein_id", "protein_length", "signalp_prediction",
    "signalp_cleavage_after", "deeptmhmm_prediction", "tm_helix_count",
    "tm_helix_coordinates", "tm_after_signal_peptide", "secretome_candidate",
    "deeptmhmm_source",
]
all_rows = []
secretome_total = 0

for fasta in sorted(proteome_dir.glob("*.all_proteins.faa")):
    sample = fasta.name.removesuffix(".all_proteins.faa")
    sequences = dict(read_fasta(fasta))
    signalp_path = signalp_dir / sample / f"{sample}_summary.signalp5"
    signalp = {}
    with signalp_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            prediction = columns[1]
            full_prediction = " ".join(columns[1:])
            cleavage_after = 0
            if prediction == "SP(Sec/SPI)":
                match = re.search(r"CS pos:\s*(\d+)-(\d+)", full_prediction)
                if match:
                    cleavage_after = int(match.group(1))
            signalp[columns[0]] = (prediction, full_prediction, cleavage_after)

    sample_rows = []
    secretome_ids = set()
    for protein_id, sequence in sequences.items():
        if protein_id not in deeptmhmm:
            raise SystemExit(f"ERROR: missing DeepTMHMM result for {protein_id}")
        if protein_id not in signalp:
            raise SystemExit(f"ERROR: missing SignalP result for {protein_id}")

        sp_class, sp_text, cleavage_after = signalp[protein_id]
        deep = deeptmhmm[protein_id]
        tm_after_sp = any(start > cleavage_after for start, _ in deep["tm_segments"])
        candidate = sp_class == "SP(Sec/SPI)" and not tm_after_sp
        if candidate:
            secretome_ids.add(protein_id)

        row = {
            "sample_id": sample,
            "protein_id": protein_id,
            "protein_length": len(sequence),
            "signalp_prediction": sp_text,
            "signalp_cleavage_after": cleavage_after or "",
            "deeptmhmm_prediction": deep["prediction"],
            "tm_helix_count": len(deep["tm_segments"]),
            "tm_helix_coordinates": ";".join(f"{start}-{end}" for start, end in deep["tm_segments"]),
            "tm_after_signal_peptide": "TRUE" if tm_after_sp else "FALSE",
            "secretome_candidate": "TRUE" if candidate else "FALSE",
            "deeptmhmm_source": deep["source_file"],
        }
        sample_rows.append(row)
        all_rows.append(row)

    sample_table = output_dir / f"{sample}.deeptmhmm_secretome.tsv"
    with sample_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_rows)

    secretome_fasta = output_dir / f"{sample}.final_secretome.faa"
    with secretome_fasta.open("w", newline="\n", encoding="utf-8") as handle:
        for protein_id, sequence in sequences.items():
            if protein_id in secretome_ids:
                handle.write(f">{protein_id}\n")
                handle.write("\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)) + "\n")
    secretome_total += len(secretome_ids)
    print(f"{sample}: {len(secretome_ids)} final secretome candidates")

master = output_dir / "deeptmhmm_secretome_master.tsv"
with master.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(all_rows)

expected = sum(1 for fasta in proteome_dir.glob("*.all_proteins.faa") for _ in read_fasta(fasta))
if len(all_rows) != expected:
    raise SystemExit(f"ERROR: expected {expected} rows but created {len(all_rows)}")

print(f"DeepTMHMM proteins integrated: {len(all_rows)}")
print(f"Final secretome candidates: {secretome_total}")
print(f"Master table: {master}")

