#!/usr/bin/env python3
"""Build corrected functional-annotation and small-protein tables."""

import argparse
import csv
import re
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--root", required=True)
parser.add_argument("--samples", required=True)
args = parser.parse_args()

root = Path(args.root)
with open(args.samples, encoding="utf-8") as handle:
    samples = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]

master_fields = [
    "sample_id", "protein_id", "protein_length", "short_protein_flag",
    "eggnog_description", "GO_terms", "KEGG_ko", "pfam_domains",
    "cazy_any_evidence", "cazy_families", "cazy_number_of_tools",
    "cazy_recommended_result", "signalp_prediction",
]
small_fields = master_fields + ["protein_sequence"]
family_pattern = re.compile(r"(?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?")


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


def nonempty(value):
    return "" if value in (None, "-") else value


normalized_dir = root / "05_normalized"
normalized_dir.mkdir(parents=True, exist_ok=True)
master_rows = []
small_rows = []

for sample in samples:
    fasta = root / "00_qc" / "fasta" / f"{sample}.faa"
    sequences = dict(read_fasta(fasta))
    annotations = {}

    for protein_id, sequence in sequences.items():
        length = len(sequence)
        annotations[protein_id] = {
            "sample_id": sample,
            "protein_id": protein_id,
            "protein_length": length,
            "short_protein_flag": "TRUE" if length < 100 else "FALSE",
            "eggnog_description": "",
            "GO_terms": "",
            "KEGG_ko": "",
            "pfam_domains": "",
            "cazy_any_evidence": "",
            "cazy_families": "",
            "cazy_number_of_tools": "",
            "cazy_recommended_result": "",
            "signalp_prediction": "",
        }

    eggnog_files = list((root / "01_eggnog" / sample).glob("*.emapper.annotations"))
    if eggnog_files:
        header = None
        with eggnog_files[0].open(encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#query"):
                    header = line[1:].rstrip("\n").split("\t")
                    continue
                if not header or line.startswith("#"):
                    continue
                record = dict(zip(header, line.rstrip("\n").split("\t")))
                query = record.get("query")
                if query in annotations:
                    annotations[query]["eggnog_description"] = nonempty(record.get("Description"))
                    annotations[query]["GO_terms"] = nonempty(record.get("GOs"))
                    annotations[query]["KEGG_ko"] = nonempty(record.get("KEGG_ko"))

    pfam_file = root / "02_pfam" / sample / f"{sample}.pfam.domtblout"
    if pfam_file.exists():
        hits = {}
        with pfam_file.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                columns = line.split()
                hits.setdefault(columns[3], set()).add(columns[0])
        for query, domains in hits.items():
            if query in annotations:
                annotations[query]["pfam_domains"] = ";".join(sorted(domains))

    overview = root / "03_dbcan" / sample / "overview.tsv"
    if overview.exists():
        with overview.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for record in reader:
                query = record.get("Gene ID", "")
                if query not in annotations:
                    continue
                tools_text = record.get("#ofTools", "0")
                try:
                    number_of_tools = int(tools_text)
                except ValueError:
                    number_of_tools = 0
                evidence_text = ";".join([
                    record.get("dbCAN_hmm", ""),
                    record.get("dbCAN_sub", ""),
                    record.get("DIAMOND", ""),
                ])
                any_families = sorted(set(family_pattern.findall(evidence_text)))
                recommended = nonempty(record.get("Recommend Results", ""))

                annotations[query]["cazy_any_evidence"] = ";".join(any_families)
                annotations[query]["cazy_number_of_tools"] = number_of_tools
                annotations[query]["cazy_recommended_result"] = recommended

                if number_of_tools >= 2:
                    annotations[query]["cazy_families"] = ";".join(any_families)

    signalp_file = root / "04_signalp" / sample / f"{sample}_summary.signalp5"
    if signalp_file.exists():
        with signalp_file.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                columns = line.rstrip("\n").split("\t")
                query = columns[0]
                if query in annotations:
                    annotations[query]["signalp_prediction"] = " ".join(columns[1:])

    sample_rows = list(annotations.values())
    sample_output = normalized_dir / f"{sample}.functional.normalized.tsv"
    with sample_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=master_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_rows)

    master_rows.extend(sample_rows)
    for row in sample_rows:
        if row["short_protein_flag"] == "TRUE":
            small_row = dict(row)
            small_row["protein_sequence"] = sequences[row["protein_id"]]
            small_rows.append(small_row)

master_dir = root / "07_master"
master_dir.mkdir(parents=True, exist_ok=True)

with (master_dir / "functional_annotation_master.tsv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=master_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(master_rows)

with (master_dir / "small_proteins_under_100aa.tsv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=small_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(small_rows)

print(f"Master proteins: {len(master_rows)}")
print(f"Proteins under 100 aa: {len(small_rows)}")
print(f"Consensus CAZymes: {sum(bool(row['cazy_families']) for row in master_rows)}")

