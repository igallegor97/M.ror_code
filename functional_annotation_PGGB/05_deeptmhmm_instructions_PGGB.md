# DeepTMHMM execution
Use the official DTU/BioLib implementation or a licensed institutional local package.
Input: `functional_annotation_C26/00_qc/MrorC26_proteins.cleaned.fa`.
Retain the full output directory, especially `predicted_topologies.3line` when produced.
Then normalize it:
```bash
python 05_parse_deeptmhmm.py   --three-line /path/to/predicted_topologies.3line   --output /Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_C26/05_deeptmhmm/deeptmhmm_normalized.tsv
```
