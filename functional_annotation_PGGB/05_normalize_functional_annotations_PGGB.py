#!/usr/bin/env python3
"""Normalize native annotation outputs into per-sample and master TSV files."""

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

fields = [
    "sample_id",
    "protein_id",
    "eggnog_description",
    "GO_terms",
    "KEGG_ko",
    "pfam_domains",
    "cazy_families",
    "signalp_prediction",
]

normalized_dir = root / "05_normalized"
normalized_dir.mkdir(parents=True, exist_ok=True)
master_rows = []

for sample in samples:
    fasta = root / "00_qc" / "fasta" / f"{sample}.faa"
    protein_ids = [
        line[1:].split()[0]
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    annotations = {
        protein_id: {
            "sample_id": sample,
            "protein_id": protein_id,
            "eggnog_description": "",
            "GO_terms": "",
            "KEGG_ko": "",
            "pfam_domains": "",
            "cazy_families": "",
            "signalp_prediction": "",
        }
        for protein_id in protein_ids
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
                    annotations[query]["eggnog_description"] = record.get("Description", "")
                    annotations[query]["GO_terms"] = record.get("GOs", "")
                    annotations[query]["KEGG_ko"] = record.get("KEGG_ko", "")

    pfam_file = root / "02_pfam" / sample / f"{sample}.pfam.domtblout"
    if pfam_file.exists():
        hits = {}
        for line in pfam_file.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("#"):
                continue
            columns = line.split()
            hits.setdefault(columns[3], set()).add(columns[0])
        for query, domains in hits.items():
            if query in annotations:
                annotations[query]["pfam_domains"] = ";".join(sorted(domains))

    dbcan_dir = root / "03_dbcan" / sample
    for candidate_name in ("overview.tsv", "overview.txt", "CAZyme_annotation.tsv"):
        candidate = dbcan_dir / candidate_name
        if not candidate.exists():
            continue
        with candidate.open(encoding="utf-8") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        for row in rows[1:]:
            if not row:
                continue
            query = row[0]
            families = set(re.findall(r"(?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?", ";".join(row)))
            if query in annotations:
                annotations[query]["cazy_families"] = ";".join(sorted(families))
        break

    signalp_dir = root / "04_signalp" / sample
    signalp_candidates = (
        list(signalp_dir.glob("*summary*"))
        + list(signalp_dir.glob("*prediction*"))
        + list(signalp_dir.glob("*.txt"))
    )
    for candidate in signalp_candidates:
        for line in candidate.read_text(errors="ignore").splitlines():
            if not line or line.startswith("#"):
                continue
            columns = line.split()
            if columns and columns[0] in annotations:
                annotations[columns[0]]["signalp_prediction"] = " ".join(columns[1:])

    sample_output = normalized_dir / f"{sample}.functional.normalized.tsv"
    with sample_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotations.values())
    master_rows.extend(annotations.values())

master_output = root / "07_master" / "functional_annotation_master.tsv"
master_output.parent.mkdir(parents=True, exist_ok=True)
with master_output.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(master_rows)

