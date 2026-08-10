#!/usr/bin/env python3
"""Parse EffectorP 3 candidate FASTA files and classify complete secretomes."""

import argparse
import csv
import re
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--secretome-dir", required=True)
parser.add_argument("--results-dir", required=True)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()

secretome_dir = Path(args.secretome_dir)
results_dir = Path(args.results_dir)
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

samples = ("MrorB3", "MrorC26", "MrorCO8", "MrorCO84", "MrorE7")
probability_pattern = re.compile(
    r"(Apoplastic|Cytoplasmic) effector probability:\s*([0-9]*\.?[0-9]+)",
    re.IGNORECASE,
)

fields = [
    "sample_id",
    "protein_id",
    "protein_length",
    "effector_candidate",
    "effectorp_primary_class",
    "effectorp_dual_localized",
    "effectorp_apoplastic_probability",
    "effectorp_cytoplasmic_probability",
    "effectorp_header",
]


def read_fasta(path):
    records = {}
    identifier = None
    header = None
    sequence = []

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if identifier in records:
                        raise ValueError(f"{path}: duplicated identifier {identifier}")
                    records[identifier] = (header, "".join(sequence))
                header = line[1:]
                identifier = header.split("|")[0].strip().split()[0]
                sequence = []
            else:
                if identifier is None:
                    raise ValueError(f"{path}: sequence before first header")
                sequence.append(line)

    if identifier is not None:
        if identifier in records:
            raise ValueError(f"{path}: duplicated identifier {identifier}")
        records[identifier] = (header, "".join(sequence))

    return records


master_rows = []

for sample in samples:
    secretome_path = secretome_dir / f"{sample}.final_secretome.faa"
    candidate_path = results_dir / sample / "EffectorCandidates.fasta"

    if not secretome_path.is_file():
        raise SystemExit(f"ERROR: missing {secretome_path}")
    if not candidate_path.is_file():
        raise SystemExit(f"ERROR: missing {candidate_path}")

    secretome = read_fasta(secretome_path)
    candidates = read_fasta(candidate_path)

    missing_ids = set(candidates) - set(secretome)
    if missing_ids:
        raise SystemExit(
            f"ERROR: {sample} has {len(missing_ids)} candidate IDs absent "
            "from the final secretome"
        )

    rows = []
    candidate_ids = set()
    primary_counts = {"Apoplastic": 0, "Cytoplasmic": 0}
    dual_count = 0

    for protein_id, (_, secretome_sequence) in secretome.items():
        row = {
            "sample_id": sample,
            "protein_id": protein_id,
            "protein_length": len(secretome_sequence),
            "effector_candidate": "FALSE",
            "effectorp_primary_class": "Non-effector",
            "effectorp_dual_localized": "FALSE",
            "effectorp_apoplastic_probability": "",
            "effectorp_cytoplasmic_probability": "",
            "effectorp_header": "",
        }

        if protein_id in candidates:
            header, candidate_sequence = candidates[protein_id]
            if candidate_sequence != secretome_sequence:
                raise SystemExit(
                    f"ERROR: candidate sequence differs from secretome for {protein_id}"
                )

            matches = probability_pattern.findall(header)
            if not matches:
                raise SystemExit(
                    f"ERROR: unrecognized EffectorP header for {protein_id}: {header}"
                )

            probabilities = {
                prediction.capitalize(): float(probability)
                for prediction, probability in matches
            }
            primary_class = matches[0][0].capitalize()
            dual = len(probabilities) == 2

            if primary_class not in primary_counts:
                raise SystemExit(
                    f"ERROR: invalid primary class for {protein_id}: {primary_class}"
                )

            if dual:
                secondary = (
                    "Cytoplasmic" if primary_class == "Apoplastic" else "Apoplastic"
                )
                if probabilities[primary_class] < probabilities[secondary]:
                    raise SystemExit(
                        f"ERROR: primary class does not have the highest "
                        f"probability for {protein_id}"
                    )

            row.update({
                "effector_candidate": "TRUE",
                "effectorp_primary_class": primary_class,
                "effectorp_dual_localized": "TRUE" if dual else "FALSE",
                "effectorp_apoplastic_probability": probabilities.get("Apoplastic", ""),
                "effectorp_cytoplasmic_probability": probabilities.get("Cytoplasmic", ""),
                "effectorp_header": header,
            })
            candidate_ids.add(protein_id)
            primary_counts[primary_class] += 1
            dual_count += int(dual)

        rows.append(row)
        master_rows.append(row)

    if len(candidate_ids) != len(candidates):
        raise SystemExit(
            f"ERROR: parsed {len(candidate_ids)} of {len(candidates)} "
            f"candidates for {sample}"
        )

    sample_output = output_dir / f"{sample}.effectorp3.tsv"
    with sample_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"{sample}: secretome={len(secretome)}, "
        f"effectors={len(candidates)}, "
        f"apoplastic_primary={primary_counts['Apoplastic']}, "
        f"cytoplasmic_primary={primary_counts['Cytoplasmic']}, "
        f"dual={dual_count}, "
        f"non_effectors={len(secretome) - len(candidates)}"
    )

master_output = output_dir / "effectorp3_master.tsv"
with master_output.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(master_rows)

print(f"Secretome proteins integrated: {len(master_rows)}")
print(
    "Effector candidates: "
    f"{sum(row['effector_candidate'] == 'TRUE' for row in master_rows)}"
)
print(f"Master table: {master_output}")

