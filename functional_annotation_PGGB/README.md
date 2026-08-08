# Functional annotation of the PacBio proteomes used in the PGGB pangenome

This directory contains the functional-annotation component of the *Moniliophthora roreri* pangenome project. It is designed for the laboratory Sun Grid Engine (SGE) cluster and follows the naming conventions used throughout the [`M.ror_code`](https://github.com/igallegor97/M.ror_code) repository.

The workflow annotates predicted proteins from the five chromosome-scale PacBio assemblies with complementary evidence from eggNOG-mapper, Pfam-A/HMMER, run_dbCAN and SignalP 5. DeepTMHMM, EffectorP 3 and selective InterProScan are run externally using FASTA files prepared by the workflow.

## Scope and relationship to PGGB

This workflow does **not** read PGGB graph files directly. GFA, VG, ODGI, partition FASTA, community VCF and SNP tables are not primary inputs. The direct input is five predicted **protein FASTA files**, one for each PacBio genome used to construct and interpret the PGGB pangenome.

The connection to PGGB is downstream:

1. PGGB defines graph communities, variants and genomic regions.
2. Population and geography analyses prioritize communities, regions and genes.
3. This workflow annotates the PacBio protein catalogues.
4. Stable identifiers connect the functional annotations to PGGB-derived candidate-gene tables.

The same annotation catalog can therefore support several PGGB analyses without repeating expensive searches. Correct mapping among GFF3 gene IDs, transcript IDs, protein FASTA IDs and PGGB candidate tables is essential.

## Directory layout

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

The layout matches the repository convention: numbered stages, descriptive names, `_PGGB` suffixes and `.sge.sh` for SGE jobs. Results and logs are written outside the source scripts according to the environment configuration.

## Required input files

### Five predicted proteomes

The biological inputs are five amino-acid FASTA files corresponding to the chromosome-scale PacBio assemblies. The workflow does not predict genes or proteins; gene prediction must already be complete.

Each FASTA must satisfy the following requirements:

- protein rather than nucleotide sequences;
- unique, non-empty record identifiers;
- the first whitespace-delimited header field is used as `protein_id`;
- permitted ID characters are letters, digits, underscore, period, colon, vertical bar, plus and minus;
- no empty sequences or unsupported residue characters.

### PacBio proteome manifest

Edit `pacbio_proteomes_manifest_PGGB.tsv`:

```text
sample_id	proteome_fasta
PacBio_1	/absolute/path/to/PacBio_1.proteins.faa
PacBio_2	/absolute/path/to/PacBio_2.proteins.faa
PacBio_3	/absolute/path/to/PacBio_3.proteins.faa
PacBio_4	/absolute/path/to/PacBio_4.proteins.faa
PacBio_5	/absolute/path/to/PacBio_5.proteins.faa
```

Replace every provisional name and path. The file must remain tab-delimited, retain its header, contain exactly five data rows, use absolute FASTA paths and contain unique sample IDs made of letters, digits, periods, underscores or hyphens.

### Optional priority-gene list

`priority_gene_ids_PGGB.txt` accepts one normalized protein ID per line. Blank lines and lines beginning with `#` are ignored. It is used to construct a selective InterProScan FASTA.

Suitable candidates include top geographic genes, genes overlapping prioritized SNPs, hotspot genes, predicted secreted proteins, CAZymes and candidate effectors.

## Software and shared databases

The workflow targets the verified laboratory HPC resources:

| Component | Version or path | Purpose |
|---|---|---|
| Python | Python 3 | FASTA validation and table normalization |
| SGE | `qsub`, arrays, `hold_jid` | scheduling and dependency control |
| eggNOG-mapper | `eggnog-mapper/2.1.13` | orthology-based annotation |
| eggNOG data | `/Storage/databases/eggnog-mapper-2.1.13` | DIAMOND proteins and annotation databases |
| HMMER | `Hmmer/3.2.1` | profile-HMM searches |
| Pfam-A | `/Storage/databases/Pfam/v35/Pfam-A.hmm` | protein families and domains |
| run_dbCAN | `/Storage/progs/run_dbCAN_4.0` | CAZyme annotation |
| dbCAN data | `/Storage/databases/dbCAN_v5.1.2` | dbCAN and dbCAN-sub HMMs |
| SignalP | `signalp/5.0b` | eukaryotic signal peptides |
| DeepTMHMM | external service | transmembrane topology |
| EffectorP | external service, version 3 | fungal effector candidates |
| InterProScan | selective external or Galaxy run | supporting signatures and ontologies |

Databases are read in place and are not duplicated in the project.

## Configuration

Edit `pipeline_environment_PGGB.env` before running the workflow.

| Variable | Meaning |
|---|---|
| `PROJECT_ROOT` | absolute path to this directory on the HPC |
| `SAMPLES_TSV` | path to the five-proteome manifest |
| `RESULTS` | output root |
| `LOGS` | SGE log directory |
| `EGGNOG_MODULE` | verified eggNOG module |
| `EGGNOG_DATA` | shared eggNOG database directory |
| `HMMER_MODULE` | verified HMMER module, including capitalization |
| `PFAM_HMM` | uncompressed and indexed Pfam-A HMM |
| `DBCAN_ENV` | run_dbCAN 4 environment |
| `DBCAN_DB` | shared dbCAN database |
| `SIGNALP_MODULE` | verified SignalP module |
| `THREADS` | default CPU count |
| `QUEUE` | optional SGE queue; blank allows automatic selection |

The `#$ -pe smp` directives determine requested SGE slots. If the cluster requires explicit `h_vmem` or `h_rt` resources, add the appropriate local SGE directives.

## Preflight validation

Run from this directory:

```bash
cd /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB
bash 00_preflight_functional_annotation_PGGB.sh
```

The preflight checks Python 3, the environment-module command, the five-row manifest, all FASTA paths, tool executables and the principal shared database files. It loads the exact verified modules and prints their versions.

The earlier observation that `emapper.py` was absent from the base shell is expected: it becomes available after loading `eggnog-mapper/2.1.13`.

Do not submit the complete workflow until the preflight ends with:

```text
OK: five samples, modules and databases verified.
```

## Running the workflow

```bash
bash run_functional_annotation_pipeline_PGGB.sh
```

The driver submits this dependency graph:

```text
proteome QC (five-task array)
        |
        +--> eggNOG (five-task array) --+
        +--> Pfam  (five-task array) ---+
        +--> dbCAN (five-task array) ---+--> integration and external FASTA preparation
        +--> SignalP (five-task array) -+
```

The printed job IDs can be monitored with `qstat`. SGE task 1 corresponds to the first manifest row, task 2 to the second, and so on.

## Methods and outputs

### Step 00 — Proteome validation and QC

Scripts:

- `00_run_proteome_qc_PGGB.sge.sh`
- `00_validate_proteomes_PGGB.py`

The validator parses each FASTA without third-party Python packages. Headers are reduced to their first field, sequences are uppercased and wrapped at 60 characters, and duplicate IDs, empty sequences and unsupported characters are rejected.

Outputs per sample:

- `results/00_qc/fasta/<sample>.faa`: normalized FASTA;
- `results/00_qc/id_maps/<sample>.id_map.tsv`: normalized ID and original header;
- `results/00_qc/stats/<sample>.qc.tsv`: protein count, amino-acid count, length statistics and invalid-residue count.

The ID map is the audit trail for reconciling annotations with original gene models.

### Step 01 — eggNOG-mapper 2.1.13

Script: `01_run_eggnog_mapper_PGGB.sge.sh`

`emapper.py` runs in protein mode with DIAMOND against the shared eggNOG database. Resume behavior permits compatible interrupted runs to reuse their intermediate results.

Principal output:

```text
results/01_eggnog/<sample>/<sample>.emapper.annotations
```

Native seed-ortholog, hits and orthology outputs may also be produced. They are retained because they contain more evidence than the compact master table.

### Step 02 — Pfam-A v35 with HMMER 3.2.1

Script: `02_run_hmmscan_pfam_PGGB.sge.sh`

`hmmscan` searches the normalized proteome against Pfam-A. The workflow uses Pfam model-specific gathering thresholds (`--cut_ga`) and produces:

- `results/02_pfam/<sample>/<sample>.pfam.domtblout`;
- `results/02_pfam/<sample>/<sample>.pfam.txt`.

The normalized result records significant model names. Use `domtblout` for domain coordinates, scores, independent E-values and repeated-domain architectures.

### Step 03 — run_dbCAN 4.0

Script: `03_run_dbcan_PGGB.sge.sh`

The installed interface uses the `CAZyme_annotation` subcommand. Because run_dbCAN 4 minor releases may expose different option names, the wrapper reads the installed help and selects supported input, output, database and CPU options. It stops rather than silently using incompatible v3 syntax.

Native results are written under:

```text
results/03_dbcan/<sample>/
```

The normalizer extracts recognized CAZy families such as GH, GT, PL, CE, AA and CBM. Consult native evidence when applying consensus rules or interpreting substrate predictions.

### Step 04 — SignalP 5.0b

Script: `04_run_signalp5_PGGB.sge.sh`

SignalP runs with eukaryotic organism mode, short output, GFF3 output, mature-protein FASTA generation, a job-local temporary directory and a controlled batch size.

Outputs are stored under:

```text
results/04_signalp/<sample>/
```

A signal-peptide prediction is not by itself a definitive secretome assignment. Proteins with additional transmembrane helices must be evaluated with DeepTMHMM.

### Step 05 — Normalization and integration

Scripts:

- `05_normalize_functional_annotations_PGGB.py`
- `05_integrate_functional_annotations_PGGB.sge.sh`
- `05_prepare_external_services_PGGB.sh`

The normalizer begins with every protein in every validated FASTA and left-joins available evidence. Unannotated proteins therefore remain in the final table.

Outputs:

- `results/05_normalized/<sample>.functional.normalized.tsv`;
- `results/07_master/functional_annotation_master.tsv`;
- service-ready FASTA files under `results/06_external/`.

The master key is `(sample_id, protein_id)` because IDs are not assumed to be globally unique across assemblies.

### Step 06 — External services

See `06_external_services_instructions_PGGB.md` for operational details.

#### DeepTMHMM

Upload each `results/06_external/<sample>.all_proteins.faa` and save all downloaded results under:

```text
results/06_external/deeptmhmm/<sample>/
```

For a conservative classical secretome, retain SignalP-positive proteins and exclude proteins with transmembrane helices after the signal-peptide region. Absence of a transmembrane helix alone does not imply secretion.

#### EffectorP 3

Upload each `<sample>.signalp_secreted.faa`, selecting fungal mode when requested. Preserve probabilities and classes, not only a binary candidate list, under:

```text
results/06_external/effectorp3/<sample>/
```

EffectorP provides candidate rankings and is not experimental proof of effector activity.

#### Selective InterProScan or Galaxy Europe

Place prioritized IDs in `priority_gene_ids_PGGB.txt` and run:

```bash
bash 05_prepare_external_services_PGGB.sh
```

Upload `results/06_external/priority_candidates.faa` to InterPro or Galaxy Europe. Request TSV and GFF3 outputs, record tool and member-database versions, and preserve the Galaxy history or workflow. Pfam has already been run locally, so additional databases provide the greatest complementary value.

InterProScan is a selective validation layer and does not block creation of the master annotation table.

## Complete output tree

```text
results/
├── 00_qc/
│   ├── fasta/
│   ├── id_maps/
│   └── stats/
├── 01_eggnog/<sample>/
├── 02_pfam/<sample>/
├── 03_dbcan/<sample>/
├── 04_signalp/<sample>/
├── 05_normalized/
│   └── <sample>.functional.normalized.tsv
├── 06_external/
│   ├── <sample>.all_proteins.faa
│   ├── <sample>.signalp_secreted.faa
│   ├── priority_candidates.faa
│   ├── deeptmhmm/
│   ├── effectorp3/
│   └── interproscan/
└── 07_master/
    └── functional_annotation_master.tsv
```

See `OUTPUT_SCHEMA_PGGB.md` for normalized column definitions.

## Connecting annotations to PGGB candidate genes

Before joining functional results to candidate-gene tables:

1. identify the GFF3 field used as the gene ID;
2. determine whether the FASTA headers contain gene, transcript or protein IDs;
3. use the QC ID maps and, if required, create a GFF3-derived gene–transcript–protein crosswalk;
4. retain sample context when combining the five annotations;
5. report unmatched and one-to-many mappings rather than silently discarding them.

This identifier-controlled join is where functional annotation becomes part of the PGGB geography analysis. Graph files are not required for the annotation searches themselves.

## Partial reruns

Individual stages can be resubmitted:

```bash
qsub -t 1-5 01_run_eggnog_mapper_PGGB.sge.sh
qsub -t 3 04_run_signalp5_PGGB.sge.sh
qsub 05_integrate_functional_annotations_PGGB.sge.sh
```

Before rerunning a stage, determine whether its output directory contains a complete compatible run. Do not combine databases or software versions without recording the difference.

## Reproducibility record

Preserve the following with the thesis analysis:

- checksums of the five input proteomes;
- the manifest and environment file;
- SGE scripts, submission date, job IDs and logs;
- software/module versions printed by the preflight;
- shared database paths and release identifiers;
- all native and normalized outputs;
- DeepTMHMM, EffectorP and InterPro/Galaxy downloads;
- external-service access dates, visible versions and selected options;
- the identifier crosswalk used to connect proteins with PGGB-derived genes.

External services can change independently of this repository, so their downloaded outputs and run metadata are necessary for reproducibility.

## Limitations

- All functional assignments are computational predictions.
- Orthology transfer can produce broad or incomplete descriptions.
- Pfam signatures do not necessarily define complete protein function.
- dbCAN family membership does not prove catalytic activity or substrate specificity.
- SignalP and DeepTMHMM describe a computational classical secretome and may miss non-classical secretion.
- EffectorP prioritizes candidates but does not demonstrate virulence.
- Annotation completeness depends on the original gene predictions.
- The compact master table does not replace the detailed native outputs.

## Minimal execution checklist

```bash
cd /Storage/data1/isabella.gallego/MAESTRIA/code/functional_annotation_PGGB

# Edit:
#   pipeline_environment_PGGB.env
#   pacbio_proteomes_manifest_PGGB.tsv

bash 00_preflight_functional_annotation_PGGB.sh
bash run_functional_annotation_pipeline_PGGB.sh
qstat
```

After completion, confirm that the master-table row count equals the sum of the five QC protein counts. Then run the external services and perform the identifier-validated join to PGGB candidate genes.
