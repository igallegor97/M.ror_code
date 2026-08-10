#!/usr/bin/env python3
"""Integrate local annotation, DeepTMHMM and EffectorP 3 into final tables."""

import argparse
import csv
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--functional-master", required=True)
parser.add_argument("--deeptmhmm-master", required=True)
parser.add_argument("--effectorp-master", required=True)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()

functional_path = Path(args.functional_master)
deep_path = Path(args.deeptmhmm_master)
effector_path = Path(args.effectorp_master)
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)


def read_table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames, rows


functional_fields, functional_rows = read_table(functional_path)
deep_fields, deep_rows = read_table(deep_path)
effector_fields, effector_rows = read_table(effector_path)

required_functional = {
    "sample_id", "protein_id", "protein_length", "short_protein_flag",
    "signalp_prediction", "cazy_families",
}
required_deep = {
    "sample_id", "protein_id", "signalp_cleavage_after",
    "deeptmhmm_prediction", "deeptmhmm_signal_length", "tm_helix_count",
    "tm_helix_coordinates", "tm_after_signal_peptide",
    "beta_barrel_prediction", "secretome_candidate", "exclusion_reason",
}
required_effector = {
    "sample_id", "protein_id", "effector_candidate",
    "effectorp_primary_class", "effectorp_dual_localized",
    "effectorp_apoplastic_probability",
    "effectorp_cytoplasmic_probability", "effectorp_header",
}

for label, required, observed in (
    ("functional", required_functional, set(functional_fields or [])),
    ("DeepTMHMM", required_deep, set(deep_fields or [])),
    ("EffectorP", required_effector, set(effector_fields or [])),
):
    missing = required - observed
    if missing:
        raise SystemExit(f"ERROR: {label} table is missing columns: {sorted(missing)}")


def key(row):
    return row["sample_id"], row["protein_id"]


records = {}
for row in functional_rows:
    row_key = key(row)
    if row_key in records:
        raise SystemExit(f"ERROR: duplicated functional key: {row_key}")

    signalp_positive = row["signalp_prediction"].startswith("SP(Sec/SPI)")
    record = dict(row)
    record.update({
        "signalp_positive": "TRUE" if signalp_positive else "FALSE",
        "deeptmhmm_tested": "FALSE",
        "signalp_cleavage_after": "",
        "deeptmhmm_prediction": "",
        "deeptmhmm_signal_length": "",
        "tm_helix_count": "",
        "tm_helix_coordinates": "",
        "tm_after_signal_peptide": "",
        "beta_barrel_prediction": "",
        "secretome_candidate": "FALSE",
        "secretome_exclusion_reason": (
            "" if signalp_positive else "SignalP_negative"
        ),
        "effectorp_tested": "FALSE",
        "effector_candidate": "",
        "effectorp_primary_class": "",
        "effectorp_dual_localized": "",
        "effectorp_apoplastic_probability": "",
        "effectorp_cytoplasmic_probability": "",
        "effectorp_header": "",
    })
    records[row_key] = record

deep_seen = set()
deep_copy = [
    "signalp_cleavage_after", "deeptmhmm_prediction",
    "deeptmhmm_signal_length", "tm_helix_count",
    "tm_helix_coordinates", "tm_after_signal_peptide",
    "beta_barrel_prediction", "secretome_candidate",
]

for row in deep_rows:
    row_key = key(row)
    if row_key in deep_seen:
        raise SystemExit(f"ERROR: duplicated DeepTMHMM key: {row_key}")
    if row_key not in records:
        raise SystemExit(f"ERROR: DeepTMHMM key absent from functional master: {row_key}")
    if records[row_key]["signalp_positive"] != "TRUE":
        raise SystemExit(f"ERROR: DeepTMHMM tested a SignalP-negative protein: {row_key}")

    deep_seen.add(row_key)
    record = records[row_key]
    record["deeptmhmm_tested"] = "TRUE"
    for column in deep_copy:
        record[column] = row[column]
    record["secretome_exclusion_reason"] = row["exclusion_reason"]

effector_seen = set()
effector_copy = [
    "effector_candidate", "effectorp_primary_class",
    "effectorp_dual_localized", "effectorp_apoplastic_probability",
    "effectorp_cytoplasmic_probability", "effectorp_header",
]

for row in effector_rows:
    row_key = key(row)
    if row_key in effector_seen:
        raise SystemExit(f"ERROR: duplicated EffectorP key: {row_key}")
    if row_key not in records:
        raise SystemExit(f"ERROR: EffectorP key absent from functional master: {row_key}")
    if records[row_key]["secretome_candidate"] != "TRUE":
        raise SystemExit(f"ERROR: EffectorP tested a non-secretome protein: {row_key}")

    effector_seen.add(row_key)
    record = records[row_key]
    record["effectorp_tested"] = "TRUE"
    for column in effector_copy:
        record[column] = row[column]

if len(records) != 104016:
    raise SystemExit(f"ERROR: expected 104016 functional proteins; found {len(records)}")
if len(deep_seen) != 8772:
    raise SystemExit(f"ERROR: expected 8772 DeepTMHMM proteins; found {len(deep_seen)}")
if len(effector_seen) != 7552:
    raise SystemExit(f"ERROR: expected 7552 EffectorP proteins; found {len(effector_seen)}")

ordered_rows = list(records.values())
secretome_rows = [row for row in ordered_rows if row["secretome_candidate"] == "TRUE"]
effector_rows_final = [row for row in ordered_rows if row["effector_candidate"] == "TRUE"]
short_rows = [row for row in ordered_rows if row["short_protein_flag"] == "TRUE"]
consensus_cazy = [row for row in ordered_rows if row["cazy_families"]]

if len(secretome_rows) != 7552:
    raise SystemExit(f"ERROR: expected 7552 secretome proteins; found {len(secretome_rows)}")
if len(effector_rows_final) != 2545:
    raise SystemExit(f"ERROR: expected 2545 effector candidates; found {len(effector_rows_final)}")
if len(short_rows) != 6109:
    raise SystemExit(f"ERROR: expected 6109 short proteins; found {len(short_rows)}")
if len(consensus_cazy) != 3553:
    raise SystemExit(f"ERROR: expected 3553 consensus CAZymes; found {len(consensus_cazy)}")

added_fields = [
    "signalp_positive", "deeptmhmm_tested", "signalp_cleavage_after",
    "deeptmhmm_prediction", "deeptmhmm_signal_length", "tm_helix_count",
    "tm_helix_coordinates", "tm_after_signal_peptide",
    "beta_barrel_prediction", "secretome_candidate",
    "secretome_exclusion_reason", "effectorp_tested",
    "effector_candidate", "effectorp_primary_class",
    "effectorp_dual_localized", "effectorp_apoplastic_probability",
    "effectorp_cytoplasmic_probability", "effectorp_header",
]
output_fields = functional_fields + added_fields


def write_table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


write_table(output_dir / "functional_annotation_master_final.tsv", ordered_rows)
write_table(output_dir / "secretome_annotations_final.tsv", secretome_rows)
write_table(output_dir / "effector_candidates_annotations_final.tsv", effector_rows_final)

summary_path = output_dir / "functional_annotation_final_summary.tsv"
summary_fields = ["sample_id", "metric", "count"]
summary_rows = []

for sample in ("MrorB3", "MrorC26", "MrorCO8", "MrorCO84", "MrorE7", "ALL"):
    subset = (
        ordered_rows if sample == "ALL"
        else [row for row in ordered_rows if row["sample_id"] == sample]
    )
    metrics = {
        "proteins": len(subset),
        "short_proteins_under_100aa": sum(row["short_protein_flag"] == "TRUE" for row in subset),
        "consensus_cazymes": sum(bool(row["cazy_families"]) for row in subset),
        "signalp_positive": sum(row["signalp_positive"] == "TRUE" for row in subset),
        "secretome_candidates": sum(row["secretome_candidate"] == "TRUE" for row in subset),
        "effector_candidates": sum(row["effector_candidate"] == "TRUE" for row in subset),
        "apoplastic_primary": sum(row["effectorp_primary_class"] == "Apoplastic" for row in subset),
        "cytoplasmic_primary": sum(row["effectorp_primary_class"] == "Cytoplasmic" for row in subset),
        "dual_localized": sum(row["effectorp_dual_localized"] == "TRUE" for row in subset),
    }
    for metric, count in metrics.items():
        summary_rows.append({"sample_id": sample, "metric": metric, "count": count})

with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

print("FINAL INTEGRATION COMPLETED")
print(f"Proteins: {len(ordered_rows)}")
print(f"Short proteins under 100 aa: {len(short_rows)}")
print(f"Consensus CAZymes: {len(consensus_cazy)}")
print(f"SignalP positive: {len(deep_seen)}")
print(f"Secretome candidates: {len(secretome_rows)}")
print(f"Effector candidates: {len(effector_rows_final)}")
print(f"Output directory: {output_dir}")

