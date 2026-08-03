#!/usr/bin/env python3
"""
generate_per_reference_vcfs.py

Generates one reference-oriented VCF for each genome represented in a PGGB
community graph.

For each sample, the script:

  1. Extracts its reference paths from the community GFA file.
  2. Validates and resolves those paths against the XG index.
  3. Associates the sample with its corresponding FASTA file.
  4. Runs vg deconstruct using all paths from that sample as reference paths.
  5. Writes one VCF and one log file per reference genome.

Required environment variables:

  GFA_FILE
      Path to the community GFA file.

  XG_FILE
      Path to the corresponding XG graph index.

  FASTA_DIR
      Directory containing one FASTA file per genome.

Optional environment variables:

  OUTPUT_DIR
      Directory where VCF and log files will be written.
      Default: per_reference_vcfs

  THREADS
      Number of threads passed to vg deconstruct.
      Default: 16

Expected graph path naming convention:

  SAMPLE#HAPLOTYPE#CONTIG

Example:

  C26#0#Group1

Usage:

    python generate_per_reference_vcfs.py

The environment variables are normally defined in the corresponding
SGE or SLURM submission script.
"""

import glob
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime

# =========================
# CONFIGURATION
# =========================

GFA_FILE = os.environ.get(
    "GFA_FILE"
)

XG_FILE = os.environ.get(
    "XG_FILE"
)

FASTA_DIR = os.environ.get(
    "FASTA_DIR"
)

OUTPUT_DIR = os.environ.get(
    "OUTPUT_DIR",
    "per_reference_vcfs",
)

THREADS_ENV = os.environ.get(
    "THREADS",
    "16",
)

# =========================
# CONSTANTS
# =========================

FASTA_PATTERNS = [
    "*.fa",
    "*.fasta",
    "*.fna",
    "*.fa.gz",
    "*.fasta.gz",
    "*.fna.gz",
]

# Messages containing any of these strings are treated as fatal vg errors.
VG_FATAL_PATTERNS = [
    "not found in graph",
    "does not exist in the graph",
    "invalid path",
    "no reference path",
    "cannot find path",
]

# =========================
# FUNCTIONS — VALIDATION
# =========================


def validate_configuration():
    """
    Validate required environment variables, paths, software, and thread count.

    Returns
    -------
    int
        Validated number of threads.
    """
    required_variables = {
        "GFA_FILE": GFA_FILE,
        "XG_FILE": XG_FILE,
        "FASTA_DIR": FASTA_DIR,
    }

    for variable_name, variable_value in required_variables.items():
        if not variable_value:
            print(
                f"[ERROR] Required environment variable "
                f"'{variable_name}' is not defined.",
                flush=True,
            )
            sys.exit(1)

    if not os.path.isfile(GFA_FILE):
        print(
            f"[ERROR] GFA_FILE does not exist or is not a file: {GFA_FILE}",
            flush=True,
        )
        sys.exit(1)

    if not os.path.isfile(XG_FILE):
        print(
            f"[ERROR] XG_FILE does not exist or is not a file: {XG_FILE}",
            flush=True,
        )
        sys.exit(1)

    if not os.path.isdir(FASTA_DIR):
        print(
            f"[ERROR] FASTA_DIR does not exist or is not a directory: "
            f"{FASTA_DIR}",
            flush=True,
        )
        sys.exit(1)

    if shutil.which("vg") is None:
        print(
            "[ERROR] vg is not available in PATH.",
            flush=True,
        )
        sys.exit(1)

    try:
        threads = int(THREADS_ENV)
    except ValueError:
        print(
            f"[ERROR] THREADS must be an integer: {THREADS_ENV}",
            flush=True,
        )
        sys.exit(1)

    if threads < 1:
        print(
            f"[ERROR] THREADS must be greater than zero: {threads}",
            flush=True,
        )
        sys.exit(1)

    return threads


# =========================
# FUNCTIONS — GFA PATHS
# =========================


def add_unique_path(paths_by_sample, sample, path_name):
    """
    Add a graph path to a sample while preserving insertion order.

    Parameters
    ----------
    paths_by_sample : dict
        Dictionary mapping sample names to path lists.

    sample : str
        Sample identifier.

    path_name : str
        Complete graph path name.
    """
    sample = sample.strip()
    path_name = path_name.strip()

    if not sample or not path_name:
        return

    paths_by_sample.setdefault(
        sample,
        [],
    )

    if path_name not in paths_by_sample[sample]:
        paths_by_sample[sample].append(
            path_name
        )


def extract_paths_from_gfa(gfa_file):
    """
    Extract sample paths from P and W records in a GFA file.

    Supported records
    -----------------
    P records
        The path name is read directly from the second column.

        Example:
            P    C26#0#Group1    ...

    W records
        The path name is reconstructed using:

            SAMPLE#HAPLOTYPE#SEQUENCE_ID

        Example:
            W    C26    0    Group1    ...

    Parameters
    ----------
    gfa_file : str
        Path to the GFA file.

    Returns
    -------
    dict
        Dictionary in the form:

            {sample_name: [path_name, ...]}
    """
    paths_by_sample = {}
    p_record_count = 0
    w_record_count = 0

    with open(
        gfa_file,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as gfa_handle:
        for line in gfa_handle:
            if not line.strip() or line.startswith("#"):
                continue

            columns = line.rstrip("\n").split("\t")
            record_type = columns[0]

            if record_type == "P":
                if len(columns) < 2:
                    continue

                path_name = columns[1]
                sample = path_name.split("#", 1)[0]

                add_unique_path(
                    paths_by_sample,
                    sample,
                    path_name,
                )

                p_record_count += 1

            elif record_type == "W":
                if len(columns) < 4:
                    continue

                sample = columns[1]
                haplotype = columns[2]
                sequence_id = columns[3]

                path_name = (
                    f"{sample}#{haplotype}#{sequence_id}"
                )

                add_unique_path(
                    paths_by_sample,
                    sample,
                    path_name,
                )

                w_record_count += 1

    print(
        f"      GFA path records: "
        f"P={p_record_count}, W={w_record_count}",
        flush=True,
    )

    return paths_by_sample


def get_paths_in_xg(xg_file):
    """
    Query the XG index and return its available graph paths.

    Parameters
    ----------
    xg_file : str
        Path to the XG index.

    Returns
    -------
    set or None
        Set of path names found in the XG index.

        Returns None when the query fails.
    """
    try:
        result = subprocess.run(
            [
                "vg",
                "paths",
                "-L",
                "-x",
                xg_file,
            ],
            capture_output=True,
            text=True,
            check=False,
        )

    except OSError as error:
        print(
            f"  [WARNING] Could not run vg paths: {error}",
            flush=True,
        )
        return None

    if result.returncode != 0:
        error_message = result.stderr.strip()

        print(
            "  [WARNING] Could not query paths from the XG index.",
            flush=True,
        )

        if error_message:
            print(
                f"            {error_message}",
                flush=True,
            )

        return None

    return {
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip()
    }


def resolve_xg_path(gfa_path, xg_paths):
    """
    Resolve a GFA path name to its exact representation in the XG index.

    Resolution order:

      1. Exact match.
      2. Match after adding a final '#0' suffix.
      3. Unique path beginning with '<GFA_PATH>#'.

    Parameters
    ----------
    gfa_path : str
        Path name extracted from the GFA file.

    xg_paths : set
        Path names present in the XG index.

    Returns
    -------
    tuple
        Tuple containing:

          - resolved path or None
          - status string

        Possible status values:
          exact
          suffix
          prefix
          ambiguous
          missing
    """
    if gfa_path in xg_paths:
        return gfa_path, "exact"

    suffixed_path = f"{gfa_path}#0"

    if suffixed_path in xg_paths:
        return suffixed_path, "suffix"

    prefix_matches = sorted(
        path
        for path in xg_paths
        if path.startswith(f"{gfa_path}#")
    )

    if len(prefix_matches) == 1:
        return prefix_matches[0], "prefix"

    if len(prefix_matches) > 1:
        return None, "ambiguous"

    return None, "missing"


def validate_paths_against_xg(paths_by_sample, xg_paths):
    """
    Resolve all GFA paths against the path names stored in the XG index.

    Ambiguous or missing paths are discarded rather than selecting one
    arbitrarily.

    Parameters
    ----------
    paths_by_sample : dict
        Paths extracted from the GFA.

    xg_paths : set
        Paths listed by vg paths.

    Returns
    -------
    dict
        Samples and paths successfully resolved against the XG index.
    """
    validated_paths = {}

    for sample, gfa_paths in sorted(
        paths_by_sample.items()
    ):
        resolved_paths = []
        missing_paths = []
        ambiguous_paths = []

        for gfa_path in gfa_paths:
            resolved_path, status = resolve_xg_path(
                gfa_path,
                xg_paths,
            )

            if resolved_path:
                if resolved_path not in resolved_paths:
                    resolved_paths.append(
                        resolved_path
                    )

            elif status == "ambiguous":
                ambiguous_paths.append(
                    gfa_path
                )

            else:
                missing_paths.append(
                    gfa_path
                )

        if missing_paths:
            print(
                f"  [WARNING] {sample}: paths not found in XG: "
                f"{missing_paths}",
                flush=True,
            )

        if ambiguous_paths:
            print(
                f"  [WARNING] {sample}: ambiguous XG path matches: "
                f"{ambiguous_paths}",
                flush=True,
            )

        if resolved_paths:
            validated_paths[sample] = resolved_paths

    return validated_paths


# =========================
# FUNCTIONS — FASTA FILES
# =========================


def list_fasta_files(fasta_dir):
    """
    Find all supported FASTA files in a directory.

    Parameters
    ----------
    fasta_dir : str
        Directory containing genome FASTA files.

    Returns
    -------
    list
        Sorted list of FASTA paths.
    """
    fasta_files = []

    for pattern in FASTA_PATTERNS:
        fasta_files.extend(
            glob.glob(
                os.path.join(
                    fasta_dir,
                    pattern,
                )
            )
        )

    return sorted(
        set(fasta_files)
    )


def remove_fasta_extensions(filename):
    """
    Remove common compressed and uncompressed FASTA extensions.

    Parameters
    ----------
    filename : str
        FASTA basename.

    Returns
    -------
    str
        Basename without its FASTA extension.
    """
    stem = filename

    if stem.lower().endswith(".gz"):
        stem = stem[:-3]

    for extension in (
        ".fasta",
        ".fna",
        ".fa",
    ):
        if stem.lower().endswith(extension):
            stem = stem[
                :-len(extension)
            ]
            break

    return stem


def find_fasta_for_sample(
    sample,
    fasta_files,
):
    """
    Find the FASTA file corresponding to a sample.

    Matching strategy:

      1. Exact normalized filename-token match.
      2. Match after removing a leading 'Mror' or 'Mror_' prefix.
      3. Delimited sample-name match.
      4. Conservative substring fallback.

    The method avoids matching CO8 to CO84.

    Parameters
    ----------
    sample : str
        Sample identifier.

    fasta_files : list
        Candidate FASTA paths.

    Returns
    -------
    tuple
        Tuple containing:

          - selected FASTA path or None
          - list of equally ranked ambiguous matches
    """
    sample_lower = sample.lower()
    ranked_matches = []

    for fasta_path in fasta_files:
        basename = os.path.basename(
            fasta_path
        )

        stem = remove_fasta_extensions(
            basename
        )

        stem_lower = stem.lower()

        without_mror = re.sub(
            r"^mror_?",
            "",
            stem_lower,
        )

        tokens = [
            token
            for token in re.split(
                r"[^a-z0-9]+",
                without_mror,
            )
            if token
        ]

        score = None

        # Highest-confidence match: sample is a complete token.
        if sample_lower in tokens:
            score = 100

        # Common annotation names such as MrorC26.groups.
        elif (
            without_mror == sample_lower
            or without_mror.startswith(
                f"{sample_lower}."
            )
            or without_mror.startswith(
                f"{sample_lower}_"
            )
            or without_mror.startswith(
                f"{sample_lower}-"
            )
        ):
            score = 95

        # Explicit alphanumeric boundary match in the original basename.
        elif re.search(
            rf"(?<![A-Za-z0-9])"
            rf"{re.escape(sample)}"
            rf"(?![A-Za-z0-9])",
            basename,
            flags=re.IGNORECASE,
        ):
            score = 90

        # Conservative fallback. Reject the match when the sample is
        # immediately followed by another alphanumeric character, which
        # prevents CO8 from matching CO84.
        else:
            sample_position = without_mror.find(
                sample_lower
            )

            if sample_position >= 0:
                end_position = (
                    sample_position
                    + len(sample_lower)
                )

                preceding_character = (
                    without_mror[sample_position - 1]
                    if sample_position > 0
                    else ""
                )

                following_character = (
                    without_mror[end_position]
                    if end_position < len(without_mror)
                    else ""
                )

                preceding_is_valid = (
                    not preceding_character
                    or not preceding_character.isalnum()
                )

                following_is_valid = (
                    not following_character
                    or not following_character.isalnum()
                )

                if preceding_is_valid and following_is_valid:
                    score = 80

        if score is not None:
            ranked_matches.append(
                (
                    score,
                    fasta_path,
                )
            )

    if not ranked_matches:
        return None, []

    highest_score = max(
        score
        for score, _ in ranked_matches
    )

    best_matches = sorted(
        fasta_path
        for score, fasta_path in ranked_matches
        if score == highest_score
    )

    if len(best_matches) == 1:
        return best_matches[0], []

    return None, best_matches


# =========================
# FUNCTIONS — VG DECONSTRUCT
# =========================


def identify_vg_messages(log_text):
    """
    Separate warning and fatal messages from vg standard error output.

    Parameters
    ----------
    log_text : str
        Text written by vg to standard error.

    Returns
    -------
    tuple
        Two lists containing warnings and fatal error messages.
    """
    warnings = []
    fatal_errors = []

    for raw_line in log_text.splitlines():
        line = raw_line.strip()

        if not line:
            continue

        lowercase_line = line.lower()

        if (
            lowercase_line.startswith("warning")
            or lowercase_line.startswith("[warning")
        ):
            warnings.append(
                line
            )

        if (
            lowercase_line.startswith("error")
            or lowercase_line.startswith("[error")
            or any(
                pattern in lowercase_line
                for pattern in VG_FATAL_PATTERNS
            )
        ):
            fatal_errors.append(
                line
            )

    return warnings, fatal_errors


def count_vcf_records(vcf_path):
    """
    Count non-header records in a plain-text VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF.

    Returns
    -------
    int
        Number of variant records.
    """
    record_count = 0

    with open(
        vcf_path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as vcf_handle:
        for line in vcf_handle:
            if line and not line.startswith("#"):
                record_count += 1

    return record_count


def vcf_has_header(vcf_path):
    """
    Check whether a VCF contains the required file-format header.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF.

    Returns
    -------
    bool
        True when the VCF header is present.
    """
    try:
        with open(
            vcf_path,
            "r",
            encoding="utf-8",
            errors="replace",
        ) as vcf_handle:
            for line in vcf_handle:
                if line.startswith(
                    "##fileformat=VCF"
                ):
                    return True

                if not line.startswith("#"):
                    break

    except OSError:
        return False

    return False


def run_vg_deconstruct(
    xg_file,
    reference_paths,
    sample_name,
    output_dir,
    log_dir,
    threads,
):
    """
    Run vg deconstruct using all graph paths from one sample as references.

    The command executed is equivalent to:

        vg deconstruct \
            -p <reference_path_1> \
            -p <reference_path_2> \
            ... \
            -a \
            -t <threads> \
            <graph.xg>

    The deprecated -e option is not used because reference-coordinate output
    is the default behavior in recent vg releases.

    Parameters
    ----------
    xg_file : str
        Path to the XG index.

    reference_paths : list
        Reference paths associated with the sample.

    sample_name : str
        Sample identifier used in output filenames.

    output_dir : str
        Directory for VCF files.

    log_dir : str
        Directory for vg log files.

    threads : int
        Number of threads.

    Returns
    -------
    tuple or None
        Tuple containing the VCF path and record count.

        Returns None when vg fails or produces an invalid VCF.
    """
    if not reference_paths:
        print(
            f"  [ERROR] No reference paths were provided for {sample_name}.",
            flush=True,
        )
        return None

    vcf_output = os.path.join(
        output_dir,
        f"{sample_name}_as_ref.vcf",
    )

    log_output = os.path.join(
        log_dir,
        f"{sample_name}_as_ref.log",
    )

    temporary_vcf = (
        f"{vcf_output}.tmp"
    )

    path_arguments = []

    for reference_path in reference_paths:
        path_arguments.extend(
            [
                "-p",
                reference_path,
            ]
        )

    command = (
        [
            "vg",
            "deconstruct",
        ]
        + path_arguments
        + [
            "-a",
            "-t",
            str(threads),
            xg_file,
        ]
    )

    print(
        f"\n  Reference : {sample_name}",
        flush=True,
    )
    print(
        f"  Paths     : {reference_paths}",
        flush=True,
    )
    print(
        f"  Command   : {' '.join(command)}",
        flush=True,
    )
    print(
        f"  Output    : {vcf_output}",
        flush=True,
    )
    print(
        f"  Log       : {log_output}",
        flush=True,
    )

    try:
        with open(
            temporary_vcf,
            "w",
            encoding="utf-8",
        ) as vcf_handle, open(
            log_output,
            "w",
            encoding="utf-8",
        ) as log_handle:
            result = subprocess.run(
                command,
                stdout=vcf_handle,
                stderr=log_handle,
                check=False,
            )

    except OSError as error:
        print(
            f"  [ERROR] Could not run vg deconstruct for "
            f"{sample_name}: {error}",
            flush=True,
        )

        if os.path.exists(temporary_vcf):
            os.remove(
                temporary_vcf
            )

        return None

    with open(
        log_output,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as log_handle:
        log_text = log_handle.read()

    warnings, fatal_errors = identify_vg_messages(
        log_text
    )

    for warning in warnings:
        print(
            f"  [WARNING] {warning}",
            flush=True,
        )

    if result.returncode != 0 or fatal_errors:
        print(
            f"  [ERROR] vg deconstruct failed for {sample_name}.",
            flush=True,
        )

        if result.returncode != 0:
            print(
                f"          Exit code: {result.returncode}",
                flush=True,
            )

        for error_message in fatal_errors:
            print(
                f"          {error_message}",
                flush=True,
            )

        print(
            f"          Full log: {log_output}",
            flush=True,
        )

        if os.path.exists(temporary_vcf):
            os.remove(
                temporary_vcf
            )

        return None

    if not os.path.isfile(temporary_vcf):
        print(
            f"  [ERROR] vg did not create a VCF for {sample_name}.",
            flush=True,
        )
        return None

    if os.path.getsize(temporary_vcf) == 0:
        print(
            f"  [ERROR] vg produced an empty VCF for {sample_name}.",
            flush=True,
        )

        os.remove(
            temporary_vcf
        )

        return None

    if not vcf_has_header(temporary_vcf):
        print(
            f"  [ERROR] Output does not contain a valid VCF header: "
            f"{temporary_vcf}",
            flush=True,
        )

        os.remove(
            temporary_vcf
        )

        return None

    variant_count = count_vcf_records(
        temporary_vcf
    )

    # Replace any previous output only after successful validation.
    os.replace(
        temporary_vcf,
        vcf_output,
    )

    if variant_count == 0:
        print(
            "  [WARNING] The generated VCF contains no variant records.",
            flush=True,
        )
    else:
        print(
            f"  Variant records: {variant_count:,}",
            flush=True,
        )

    return vcf_output, variant_count


# =========================
# MAIN WORKFLOW
# =========================


def main():
    """Run the complete per-reference VCF generation workflow."""
    threads = validate_configuration()

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    log_dir = os.path.join(
        OUTPUT_DIR,
        "logs",
    )

    os.makedirs(
        log_dir,
        exist_ok=True,
    )

    start_time = datetime.now()

    print("=" * 65, flush=True)
    print("PER-REFERENCE VCF GENERATION", flush=True)
    print(
        f"Start       : {start_time.strftime('%Y-%m-%d %H:%M:%S')}",
        flush=True,
    )
    print(f"GFA_FILE    : {GFA_FILE}", flush=True)
    print(f"XG_FILE     : {XG_FILE}", flush=True)
    print(f"FASTA_DIR   : {FASTA_DIR}", flush=True)
    print(f"OUTPUT_DIR  : {OUTPUT_DIR}", flush=True)
    print(f"THREADS     : {threads}", flush=True)
    print("=" * 65, flush=True)

    # ---------------------------------------------------------
    # 1. Extract and validate graph paths
    # ---------------------------------------------------------

    print(
        "\n[1/3] Reading graph paths from the GFA file...",
        flush=True,
    )

    paths_by_sample = extract_paths_from_gfa(
        GFA_FILE
    )

    if not paths_by_sample:
        print(
            "[ERROR] No P or W path records were found in the GFA file.",
            flush=True,
        )
        sys.exit(1)

    print(
        f"      Samples found: {len(paths_by_sample)}",
        flush=True,
    )

    for sample, paths in sorted(
        paths_by_sample.items()
    ):
        print(
            f"        {sample}: {len(paths)} path(s)",
            flush=True,
        )

        for path_name in paths:
            print(
                f"          {path_name}",
                flush=True,
            )

    print(
        "\n      Validating paths against the XG index...",
        flush=True,
    )

    xg_paths = get_paths_in_xg(
        XG_FILE
    )

    if xg_paths is not None:
        print(
            f"      Paths found in XG: {len(xg_paths)}",
            flush=True,
        )

        example_paths = sorted(
            xg_paths
        )[:3]

        if example_paths:
            print(
                f"      XG path examples: {example_paths}",
                flush=True,
            )

        paths_by_sample = validate_paths_against_xg(
            paths_by_sample,
            xg_paths,
        )

        if not paths_by_sample:
            print(
                "[ERROR] No GFA paths could be resolved against the XG index.",
                flush=True,
            )
            sys.exit(1)

        print(
            f"      Samples with valid XG paths: "
            f"{len(paths_by_sample)}",
            flush=True,
        )

        for sample, paths in sorted(
            paths_by_sample.items()
        ):
            print(
                f"        {sample}: {len(paths)} valid path(s)",
                flush=True,
            )

    else:
        print(
            "      [WARNING] XG path validation was unavailable. "
            "The original GFA path names will be used.",
            flush=True,
        )

    # ---------------------------------------------------------
    # 2. Associate samples with FASTA files
    # ---------------------------------------------------------

    print(
        "\n[2/3] Associating FASTA files with samples...",
        flush=True,
    )

    fasta_files = list_fasta_files(
        FASTA_DIR
    )

    if not fasta_files:
        print(
            f"[ERROR] No supported FASTA files were found in: {FASTA_DIR}",
            flush=True,
        )
        print(
            "        Supported extensions: "
            ".fa, .fasta, .fna, .fa.gz, .fasta.gz, .fna.gz",
            flush=True,
        )
        sys.exit(1)

    print(
        f"      FASTA files available: {len(fasta_files)}",
        flush=True,
    )

    for fasta_path in fasta_files:
        print(
            f"        {os.path.basename(fasta_path)}",
            flush=True,
        )

    sample_fastas = {}
    samples_without_fasta = []
    samples_with_ambiguous_fasta = []

    for sample in sorted(
        paths_by_sample
    ):
        fasta_path, ambiguous_matches = find_fasta_for_sample(
            sample,
            fasta_files,
        )

        if fasta_path:
            sample_fastas[sample] = fasta_path

            print(
                f"        {sample} → {os.path.basename(fasta_path)}",
                flush=True,
            )

        elif ambiguous_matches:
            print(
                f"  [WARNING] Multiple equally ranked FASTA files "
                f"matched sample '{sample}':",
                flush=True,
            )

            for candidate in ambiguous_matches:
                print(
                    f"            {os.path.basename(candidate)}",
                    flush=True,
                )

            samples_with_ambiguous_fasta.append(
                sample
            )

        else:
            print(
                f"  [WARNING] No FASTA file found for sample: {sample}",
                flush=True,
            )

            samples_without_fasta.append(
                sample
            )

    if not sample_fastas:
        print(
            "[ERROR] No samples could be associated with FASTA files.",
            flush=True,
        )
        sys.exit(1)

    # ---------------------------------------------------------
    # 3. Run vg deconstruct for each reference
    # ---------------------------------------------------------

    print(
        "\n[3/3] Running vg deconstruct for each reference genome...",
        flush=True,
    )
    print(
        "-" * 65,
        flush=True,
    )

    results = {}
    failed_samples = []

    for sample in sorted(
        paths_by_sample
    ):
        if sample not in sample_fastas:
            print(
                f"\n  [SKIP] {sample}: no unambiguous FASTA association.",
                flush=True,
            )

            failed_samples.append(
                sample
            )

            continue

        reference_paths = paths_by_sample[
            sample
        ]

        result = run_vg_deconstruct(
            xg_file=XG_FILE,
            reference_paths=reference_paths,
            sample_name=sample,
            output_dir=OUTPUT_DIR,
            log_dir=log_dir,
            threads=threads,
        )

        if result is None:
            failed_samples.append(
                sample
            )

        else:
            vcf_path, variant_count = result

            results[sample] = {
                "vcf": vcf_path,
                "variant_count": variant_count,
                "fasta": sample_fastas[sample],
                "path_count": len(reference_paths),
            }

    # ---------------------------------------------------------
    # Final report
    # ---------------------------------------------------------

    end_time = datetime.now()
    elapsed_time = end_time - start_time

    print("\n" + "=" * 65, flush=True)
    print("PER-REFERENCE VCF SUMMARY", flush=True)

    if results:
        print(
            "\n"
            f"{'Reference':<14} "
            f"{'Paths':>7} "
            f"{'Variants':>14}  "
            f"VCF",
            flush=True,
        )

        print(
            "-" * 65,
            flush=True,
        )

        for sample, information in sorted(
            results.items()
        ):
            print(
                f"{sample:<14} "
                f"{information['path_count']:>7} "
                f"{information['variant_count']:>14,}  "
                f"{os.path.basename(information['vcf'])}",
                flush=True,
            )

    if samples_without_fasta:
        print(
            "\nSamples without FASTA files: "
            + ", ".join(samples_without_fasta),
            flush=True,
        )

    if samples_with_ambiguous_fasta:
        print(
            "Samples with ambiguous FASTA matches: "
            + ", ".join(samples_with_ambiguous_fasta),
            flush=True,
        )

    if failed_samples:
        print(
            "Failed or skipped samples: "
            + ", ".join(sorted(set(failed_samples))),
            flush=True,
        )

    print(
        f"\nSamples detected    : {len(paths_by_sample)}",
        flush=True,
    )
    print(
        f"VCFs generated      : {len(results)}",
        flush=True,
    )
    print(
        f"VCFs failed/skipped : {len(set(failed_samples))}",
        flush=True,
    )
    print(
        f"Results saved to    : {OUTPUT_DIR}",
        flush=True,
    )
    print(
        f"Logs saved to       : {log_dir}",
        flush=True,
    )
    print(
        f"Total runtime       : {elapsed_time}",
        flush=True,
    )
    print(
        f"End                 : "
        f"{end_time.strftime('%Y-%m-%d %H:%M:%S')}",
        flush=True,
    )
    print("=" * 65, flush=True)

    print(
        "\nNext step: classify the generated VCF files with the "
        "appropriate per-reference variant-classification script.",
        flush=True,
    )

    if failed_samples and not results:
        sys.exit(1)


if __name__ == "__main__":
    main()
