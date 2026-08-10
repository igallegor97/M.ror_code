# External functional-annotation services

## DeepTMHMM

Upload each `results/06_external/<sample>.all_proteins.faa` file to the official
DeepTMHMM service. Download the complete prediction and topology files and save
them under:

```text
results/06_external/deeptmhmm/<sample>/
```

For a conservative classical secretome, combine a positive SignalP prediction
with the absence of additional transmembrane helices after the signal-peptide
region. A protein without a transmembrane helix is not automatically secreted.

## EffectorP 3

Upload each `<sample>.signalp_secreted.faa` and select the fungal protein mode
when requested. Preserve probabilities and prediction classes under:

```text
results/06_external/effectorp3/<sample>/
```

EffectorP results are computational candidates, not experimental evidence.

## Selective InterProScan or Galaxy Europe

1. Add one prioritized protein ID per line to `priority_gene_ids_PGGB.txt`.
2. Run `bash 05_prepare_external_services_PGGB.sh`.
3. Upload `results/06_external/priority_candidates.faa`.
4. Request TSV and GFF3 outputs and GO/pathway mappings when available.
5. Record the date, tool version, selected databases and parameters.
6. Save results under `results/06_external/interproscan/`.

Pfam has already been searched locally. Selecting complementary InterPro member
databases provides the greatest additional value. InterProScan is a selective
validation layer and is not required for the local master table.

