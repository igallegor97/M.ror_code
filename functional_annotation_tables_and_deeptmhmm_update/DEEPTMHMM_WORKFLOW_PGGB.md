# DeepTMHMM and secretome workflow

## 1. Regenerate corrected local tables

Replace `05_normalize_functional_annotations_PGGB.py` in `PROJECT_ROOT`, then:

```bash
cd "$PROJECT_ROOT"
source pipeline_config_PGGB.sh

python3 05_normalize_functional_annotations_PGGB.py \
  --root "$RESULTS" \
  --samples "$SAMPLES_TSV"
```

Expected summary:

```text
Master proteins: 104016
Proteins under 100 aa: 6109
Consensus CAZymes: 3553
```

## 2. Prepare DeepTMHMM batches

Copy `06_prepare_deeptmhmm_inputs_PGGB.py` into `PROJECT_ROOT`, then:

```bash
python3 06_prepare_deeptmhmm_inputs_PGGB.py \
  --input-dir "$RESULTS/06_external" \
  --output-dir "$RESULTS/06_external/deeptmhmm_inputs" \
  --batch-size 2000
```

The script stops if a sequence contains characters outside the 20-amino-acid
alphabet accepted by the official service.

## 3. Submit batches

Open:

https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/

Upload each batch listed in:

```text
$RESULTS/06_external/deeptmhmm_inputs/deeptmhmm_batch_manifest.tsv
```

Multiple-sequence jobs do not produce individual plots. Download the complete
result archive for every batch. The required file is normally named:

```text
predicted_topologies.3line
```

Rename it according to the manifest's `expected_output` column and store it under:

```text
$RESULTS/06_external/deeptmhmm_results/<sample>/
```

Do not concatenate files manually and do not alter protein identifiers.

## 4. Verify downloaded coverage

```bash
find "$RESULTS/06_external/deeptmhmm_results" \
  -type f -name '*.3line' | sort
```

There should be one output for every manifest row.

## 5. Parse and construct the final secretome

Copy `06_parse_deeptmhmm_and_build_secretome_PGGB.py` into `PROJECT_ROOT`, then:

```bash
python3 06_parse_deeptmhmm_and_build_secretome_PGGB.py \
  --results "$RESULTS/06_external/deeptmhmm_results" \
  --signalp-dir "$RESULTS/04_signalp" \
  --proteome-dir "$RESULTS/06_external" \
  --output-dir "$RESULTS/06_external/secretome_final"
```

The parser requires one DeepTMHMM result and one SignalP result for every one of
the 104,016 proteins. A secretome candidate must be SignalP-positive and have no
DeepTMHMM transmembrane segment beginning after the SignalP cleavage position.

The final EffectorP 3 inputs are:

```text
$RESULTS/06_external/secretome_final/<sample>.final_secretome.faa
```

