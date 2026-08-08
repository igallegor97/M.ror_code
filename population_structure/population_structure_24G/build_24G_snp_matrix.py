#!/usr/bin/env python3
"""
build_24G_snp_matrix.py

Builds a binary SNP matrix from selected PGGB community VCFs.

Profiles:
  all_11         - all chromosome-associated communities
  conservative_9 - excludes community.0 and community.1

The reference omitted by vg deconstruct is reconstructed as genotype 0
separately for every community.
"""

from __future__ import annotations
import argparse, gzip, re, sys
from collections import Counter
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd

def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument("--profile",required=True,choices=["all_11","conservative_9"])
    p.add_argument("--base-dir",required=True)
    p.add_argument("--reference-summary",required=True)
    p.add_argument("--manifest",required=True)
    p.add_argument("--name-map",required=True)
    p.add_argument("--metadata",required=True)
    p.add_argument("--output-dir",required=True)
    p.add_argument("--max-missing",type=float,default=0.20)
    return p.parse_args()

def natural_key(v):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)",v)]

def open_vcf(path):
    return gzip.open(path,"rt") if str(path).endswith(".gz") else open(path)

def gt_binary(field,gt_idx):
    parts=field.split(":")
    gt=parts[gt_idx] if gt_idx < len(parts) else "."
    if gt in {"",".","./.",".|."}: return np.nan
    alleles=re.split(r"[/|]",gt)
    if any(a in {"","."} for a in alleles): return np.nan
    try: nums=[int(a) for a in alleles]
    except ValueError: return np.nan
    return 1.0 if any(a>0 for a in nums) else 0.0

def load_map(path):
    df=pd.read_csv(path,sep="\t",dtype=str)
    if set(["vcf_sample","sample_id"])-set(df.columns):
        raise ValueError("Invalid sample name map.")
    if df.vcf_sample.duplicated().any() or df.sample_id.duplicated().any():
        raise ValueError("Duplicated names in sample map.")
    return dict(zip(df.vcf_sample,df.sample_id))

def load_refs(path):
    df=pd.read_csv(path,sep="\t",dtype=str)
    need={"community","reference_prefix","status"}
    if need-set(df.columns): raise ValueError("Invalid reference summary.")
    df=df[df.status=="SUCCESS"]
    return dict(zip(df.community,df.reference_prefix))

def load_communities(path,profile):
    df=pd.read_csv(path,sep="\t",dtype=str)
    col="analysis_all" if profile=="all_11" else "analysis_conservative"
    return sorted(df.loc[df[col].str.upper()=="YES","community"].tolist(),key=natural_key)

def header_samples(path):
    with open_vcf(path) as h:
        for line in h:
            if line.startswith("#CHROM"):
                return line.rstrip().split("\t")[9:]
    raise ValueError(f"No #CHROM line in {path}")

def main():
    a=parse_args()
    start=datetime.now()
    base=Path(a.base_dir)
    out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    mapping=load_map(a.name_map)
    refs=load_refs(a.reference_summary)
    comms=load_communities(a.manifest,a.profile)
    meta=pd.read_csv(a.metadata,sep="\t",dtype=str)
    samples=meta.sample_id.tolist()
    if len(samples)!=24 or len(set(samples))!=24:
        raise ValueError("Metadata must contain 24 unique samples.")

    rows=[]; feat=[]; qc=[]
    miss=Counter(); totals=Counter(); feature_no=0

    print("="*70,flush=True)
    print(f"PROFILE: {a.profile}",flush=True)
    print(f"COMMUNITIES: {','.join(comms)}",flush=True)
    print("="*70,flush=True)

    for comm in comms:
        vcf=base/f"all_genomes_pansn.fasta.bf3285f.{comm}"/"variants.vcf.gz"
        if not vcf.is_file(): raise FileNotFoundError(vcf)
        if comm not in refs: raise ValueError(f"No reference for {comm}")
        raw_ref=refs[comm]
        if raw_ref not in mapping: raise ValueError(f"Reference not in map: {raw_ref}")
        ref_sample=mapping[raw_ref]

        raw_samples=header_samples(vcf)
        unknown=[s for s in raw_samples if s not in mapping]
        if unknown: raise ValueError(f"Unknown samples in {comm}: {unknown}")
        norm=[mapping[s] for s in raw_samples]
        reconstructed = set(norm) | {ref_sample}

        # Some chromosome-associated communities do not contain all 24
        # genomes. This is expected for community.0, where B3 Group10 was
        # partitioned into community.1. Samples absent from the graph must
        # remain NA rather than being reconstructed as reference genotypes.
        extra_samples = sorted(reconstructed - set(samples))
        absent_samples = sorted(set(samples) - reconstructed)

        if extra_samples:
            raise ValueError(
                f"{comm}: graph contains samples absent from metadata: "
                f"{extra_samples}"
            )

        if absent_samples:
            print(
                f"[WARNING] {comm}: samples structurally absent from graph "
                f"and retained as NA: {absent_samples}",
                flush=True,
            )

        raw_records=0; retained_here=0
        with open_vcf(vcf) as h:
            for line in h:
                if line.startswith("#"): continue
                raw_records+=1
                f=line.rstrip().split("\t")
                if len(f)<10: continue
                chrom,pos,ref,alt=f[0],f[1],f[3],f[4]
                if "," in alt or len(ref)!=1 or len(alt)!=1: continue
                if ref in {".","*"} or alt in {".","*"}: continue
                fmt=f[8].split(":")
                if "GT" not in fmt: continue
                gt_idx=fmt.index("GT")

                vals={s:np.nan for s in samples}
                vals[ref_sample]=0.0
                for raw,field in zip(raw_samples,f[9:]):
                    vals[mapping[raw]]=gt_binary(field,gt_idx)
                vector=[vals[s] for s in samples]
                observed=[x for x in vector if not np.isnan(x)]
                if len(observed)<2 or len(set(observed))<2: continue

                feature_no+=1
                fid=f"{comm}|{chrom}|{pos}|{ref}|{alt}|{feature_no}"
                rows.append([fid]+vector)
                feat.append({
                    "feature_id":fid,"community":comm,"vcf_file":str(vcf),
                    "reference_sample":ref_sample,"chrom":chrom,"pos":int(pos),
                    "ref":ref,"alt":alt
                })
                retained_here+=1
                for s,x in zip(samples,vector):
                    totals[s]+=1
                    if np.isnan(x): miss[s]+=1

        qc.append({
            "community":comm,"vcf_file":str(vcf),
            "raw_vcf_samples":",".join(raw_samples),
            "normalized_vcf_samples":",".join(norm),
            "reference_sample_added": ref_sample,
            "samples_structurally_absent": ",".join(absent_samples),
            "n_samples_structurally_absent": len(absent_samples),
            "samples_reconstructed_total": len(reconstructed),
            "vcf_records": raw_records,
            "polymorphic_biallelic_snps": retained_here
        })
        print(f"{comm}: records={raw_records:,}; SNPs={retained_here:,}; ref={ref_sample}",flush=True)

    raw=pd.DataFrame(rows,columns=["feature_id"]+samples).set_index("feature_id")
    fm=pd.DataFrame(feat).set_index("feature_id")
    if raw.empty: raise RuntimeError("No SNPs retained.")

    missing_fraction=raw.isna().mean(axis=1)
    filt=raw.loc[missing_fraction<=a.max_missing].copy()
    fm=fm.loc[filt.index].copy()
    fm["missing_fraction"]=missing_fraction.loc[filt.index]

    means=filt.mean(axis=1)
    filt=filt.loc[means.notna()].copy()
    fm=fm.loc[filt.index].copy()
    means=means.loc[filt.index]
    imp=filt.T.fillna(means).T

    variances=imp.var(axis=1)
    imp=imp.loc[variances>0]
    filt=filt.loc[imp.index]
    fm=fm.loc[imp.index]
    fm["imputed_variance"]=variances.loc[imp.index]

    prefix=f"PGGB24_{a.profile}"
    raw.to_csv(out/f"{prefix}_snp_matrix_raw.tsv",sep="\t")
    filt.to_csv(out/f"{prefix}_snp_matrix_filtered.tsv",sep="\t")
    imp.to_csv(out/f"{prefix}_snp_matrix_imputed.tsv",sep="\t")
    fm.to_csv(out/f"{prefix}_snp_feature_metadata.tsv",sep="\t")
    pd.DataFrame(qc).to_csv(out/f"{prefix}_vcf_qc.tsv",sep="\t",index=False)

    sq=[]
    for s in samples:
        total=int(totals[s]); missing=int(miss[s])
        sq.append({
            "sample_id":s,
            "total_features_before_filter":total,
            "missing_genotypes_before_filter":missing,
            "missing_fraction_before_filter":round(missing/total,6) if total else 0
        })
    pd.DataFrame(sq).to_csv(out/f"{prefix}_sample_missingness.tsv",sep="\t",index=False)

    contribution=fm.groupby("community").size().rename("retained_features").reset_index()
    contribution.to_csv(out/f"{prefix}_community_contribution.tsv",sep="\t",index=False)

    summary=pd.DataFrame([
        ["profile",a.profile],
        ["communities",len(comms)],
        ["samples",len(samples)],
        ["raw_polymorphic_features",raw.shape[0]],
        ["features_retained_final",imp.shape[0]],
        ["max_missing",a.max_missing]
    ],columns=["metric","value"])
    summary.to_csv(out/f"{prefix}_matrix_summary.tsv",sep="\t",index=False)

    print("="*70,flush=True)
    print(f"FINAL FEATURES: {imp.shape[0]:,}",flush=True)
    print(f"SAMPLES: {imp.shape[1]}",flush=True)
    print(f"RUNTIME: {datetime.now()-start}",flush=True)
    print("="*70,flush=True)

if __name__=="__main__":
    try: main()
    except Exception as e:
        print(f"[ERROR] {e}",file=sys.stderr,flush=True)
        raise
