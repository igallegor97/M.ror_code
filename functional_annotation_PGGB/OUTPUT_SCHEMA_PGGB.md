# Normalized functional-annotation schema

The stable master key is `(sample_id, protein_id)`.

| Column | Description |
|---|---|
| `sample_id` | PacBio assembly identifier from the manifest |
| `protein_id` | First field of the normalized FASTA header |
| `eggnog_description` | eggNOG-mapper functional description |
| `GO_terms` | Gene Ontology terms transferred by eggNOG-mapper |
| `KEGG_ko` | KEGG orthology assignments from eggNOG-mapper |
| `pfam_domains` | Significant Pfam-A models detected with `--cut_ga` |
| `cazy_families` | CAZy families extracted from run_dbCAN output |
| `signalp_prediction` | SignalP 5 prediction record |

Native tool outputs remain authoritative for coordinates, scores, probabilities,
thresholds and detailed evidence.

