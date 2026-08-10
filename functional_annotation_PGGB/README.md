# Functional annotation of the PacBio proteomes used in the PGGB pangenome

This directory contains the functional-annotation workflow for the five
chromosome-scale PacBio assemblies of *Moniliophthora roreri*. It follows the
organization and English-language shell-script style used throughout the
`M.ror_code` repository: numbered stages, `_PGGB` suffixes, SGE scripts ending
in `.sge.sh`, structured header blocks and explicit execution reports.

## Scope

This workflow does not analyze PGGB GFA, VG, ODGI, VCF, community or SNP files
directly. Its direct inputs are the five predicted protein FASTA files used in
the PacBio component of the pangenome project.

The connection to PGGB is downstream:

```text
PGGB communities and variants
             |
population and geography analyses
             |
prioritized genes and genomic regions
             |
protein identifier mapping
             |
functional annotations generated here
```

The annotation catalog can therefore be reused across PGGB analyses. Gene,
transcript and protein identifiers must be reconciled before joining the master
annotation table to candidate-gene tables.

## Input proteomes

The configured manifest contains the five protein FASTA files:

| Sample | Protein FASTA |
|---|---|
| MrorB3 | `/Storage/data1/isabella.gallego/MAESTRIA/results_pacbio/MrorB3.groups_proteins.fa` |
| MrorC26 | `/Storage/data1/isabella.gallego/MAESTRIA/results_pacbio/MrorC26.groups_proteins.fa` |
| MrorCO8 | `/Storage/data1/isabella.gallego/MAESTRIA/results_pacbio/MrorCO8.groups_proteins.fa` |
| MrorCO84 | `/Storage/data1/isabella.gallego/MAESTRIA/results_pacbio/MrorCO84.groups_proteins.fa` |
| MrorE7 | `/Storage/data1/isabella.gallego/MAESTRIA/results_pacbio/MrorE7.groups_proteins.fa` |

The workflow does not perform gene prediction. Each FASTA must contain amino-acid
sequences with unique, non-empty identifiers. The first whitespace-delimited
field of each header becomes `protein_id`.

## Directory contents

```text
functional_annotation_PGGB/
├── 00_preflight_functional_annotation_PGGB.sh
├── 00_run_proteome_qc_PGGB.sge.sh
├── 00_validate_proteomes_PGGB.py
├── 01_run_eggnog_mapper_PGGB.sge.sh
├── 02_run_hmmscan_pfam_PGGB.sge.sh
├── 03_run_dbcan_PGGB.sge.sh
├── 04_run_signalp5_PGGB.sge.sh
├── 05_integrate_functional_annotations_PGGB.sge.sh
├── 05_normalize_functional_annotations_PGGB.py
├── 05_prepare_external_services_PGGB.sh
├── 06_external_services_instructions_PGGB.md
├── OUTPUT_SCHEMA_PGGB.md
├── pacbio_proteomes_manifest_PGGB.tsv
├── pipeline_config_PGGB.sh
├── pipeline_environment_PGGB.env
├── priority_gene_ids_PGGB.txt
├── run_functional_annotation_pipeline_PGGB.sh
└── README.md
```

## Configured HPC environment

| Resource | Configuration |
|---|---|
| Project directory | `/Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB` |
| Result directory | `/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_PGGB_results` |
| Log directory | `$PROJECT_ROOT/logs` |
| SGE queue | `all.q` |
| Default CPUs | 8 |
| Memory target | 32 GB |
| Runtime target | 48 hours |
| eggNOG-mapper | `eggnog-mapper/2.1.13` |
| eggNOG database | `/Storage/databases/eggnog-mapper-2.1.13` |
| HMMER | `Hmmer/3.2.1` |
| Pfam-A | `/Storage/databases/Pfam/v35/Pfam-A.hmm` |
| run_dbCAN | `/Storage/progs/run_dbCAN_4.0` |
| dbCAN database | `/Storage/databases/dbCAN_v5.1.2` |
| SignalP | `signalp/5.0b` |

The shared databases are read in place and are not copied to the project.

## Preflight validation

Run:

```bash
cd /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB
bash 00_preflight_functional_annotation_PGGB.sh
```

The preflight validates:

- Python 3 and environment modules;
- the manifest header and exactly five unique samples;
- readability of all five protein FASTA files;
- eggNOG-mapper and its principal shared databases;
- HMMER, Pfam-A and the four pressed HMM indexes;
- run_dbCAN 4 and the dbCAN HMM files;
- SignalP 5.0b.

`emapper.py` is expected to be absent before loading its module.

## Complete execution

After a successful preflight:

```bash
bash run_functional_annotation_pipeline_PGGB.sh
```

The driver submits:

```text
proteome QC (array 1-5)
       |
       +-- eggNOG (array 1-5) --+
       +-- Pfam  (array 1-5) ---+
       +-- dbCAN (array 1-5) ---+-- integration
       +-- SignalP (array 1-5) -+
```

Monitor the jobs with:

```bash
qstat
```

All job scripts are submitted through absolute paths. `PROJECT_ROOT` is exported
with `qsub -V`, so compute nodes do not depend on the SGE spool directory or the
submission working directory. The driver also converts a terse array ID such as
`25581.1-5:1` to the base job ID `25581` before using `-hold_jid`.

## Step 00: proteome QC

The validator:

- reduces each header to its first field;
- converts sequences to uppercase;
- wraps normalized FASTA sequences at 60 characters;
- rejects duplicate IDs, empty sequences and unsupported characters;
- preserves the original header in an ID map;
- reports protein counts and length statistics.

Outputs:

```text
00_qc/fasta/<sample>.faa
00_qc/id_maps/<sample>.id_map.tsv
00_qc/stats/<sample>.qc.tsv
```

## Step 01: eggNOG-mapper

eggNOG-mapper 2.1.13 runs in protein mode with DIAMOND and the shared database.
It supplies orthology-based descriptions, GO terms and KEGG orthology entries.
Native results are stored under:

```text
01_eggnog/<sample>/
```

The workflow uses `--resume` to allow compatible interrupted runs to continue.

## Step 02: Pfam-A

`hmmscan` searches every protein against Pfam-A v35 using model-specific
gathering thresholds (`--cut_ga`). Outputs are:

```text
02_pfam/<sample>/<sample>.pfam.domtblout
02_pfam/<sample>/<sample>.pfam.txt
```

The domain table retains coordinates, E-values and scores and should be used for
detailed domain-architecture interpretation.

## Step 03: run_dbCAN

run_dbCAN 4.0 annotates CAZymes in protein mode. The wrapper inspects the
installed `CAZyme_annotation --help` output and selects the supported v4 input,
output, database and CPU options rather than assuming obsolete v3 syntax.

Outputs are stored under:

```text
03_dbcan/<sample>/
```

CAZy-family membership does not by itself prove enzymatic activity or substrate
specificity.

## Step 04: SignalP

SignalP 5.0b runs in eukaryotic mode and requests short tabular, GFF3 and mature
protein outputs. Results are stored under:

```text
04_signalp/<sample>/
```

A positive SignalP result is evidence for a classical signal peptide, not a
complete secretome classification. DeepTMHMM is used to evaluate additional
transmembrane helices.

## Step 05: normalization and integration

The integration stage begins with every protein in the validated FASTA files and
left-joins the available evidence. Proteins without annotations remain present.

Outputs:

```text
05_normalized/<sample>.functional.normalized.tsv
07_master/functional_annotation_master.tsv
```

The stable master key is `(sample_id, protein_id)`.

## Step 06: external tools

The pipeline prepares:

```text
06_external/<sample>.all_proteins.faa
06_external/<sample>.signalp_secreted.faa
06_external/priority_candidates.faa
```

- Submit all-protein FASTA files to DeepTMHMM.
- Submit SignalP-secreted FASTA files to EffectorP 3.
- Add prioritized protein IDs to `priority_gene_ids_PGGB.txt` and regenerate
  `priority_candidates.faa` for selective InterProScan or Galaxy Europe.

See `06_external_services_instructions_PGGB.md`.

## Output tree

```text
/Storage/data1/isabella.gallego/MAESTRIA/functional_annotation_PGGB_results/
├── 00_qc/
├── 01_eggnog/
├── 02_pfam/
├── 03_dbcan/
├── 04_signalp/
├── 05_normalized/
├── 06_external/
└── 07_master/
```

## Joining annotations to PGGB candidate genes

Before joining these results to community-geography outputs:

1. determine whether FASTA headers represent genes, transcripts or proteins;
2. identify the corresponding fields in the reference GFF3;
3. use the QC ID maps and, when necessary, a gene–transcript–protein crosswalk;
4. retain `sample_id` when combining the five proteomes;
5. report unmatched and one-to-many mappings.

## Reproducibility record

Preserve:

- input FASTA checksums;
- this manifest and environment file;
- preflight output and software versions;
- SGE job IDs and logs;
- database paths and releases;
- all native and normalized outputs;
- dates, versions and parameters for external services;
- the identifier crosswalk used to connect annotations to PGGB genes.

## Minimal checklist

```bash
cd /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB
bash 00_preflight_functional_annotation_PGGB.sh
bash run_functional_annotation_pipeline_PGGB.sh
qstat
```

After completion, verify that the master-table row count equals the sum of the
five QC protein counts before integrating the annotations with PGGB-derived
candidate genes.

