#!/usr/bin/env python3
"""
annotate_community_snps.py

Maps retained SNPs from the final feature metadata to genes and genomic
regions in a reference GFF3. Supports flexible column names for chromosome
and position.

Outputs:
  community_snp_gene_hits.tsv
  community_region_summary.tsv
  community_gene_summary.tsv
"""

from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import pandas as pd


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument("--features",required=True)
    p.add_argument("--gff3",required=True)
    p.add_argument("--output-dir",required=True)
    p.add_argument("--reference-prefix",default="")
    return p.parse_args()


def pick_col(df,cands):
    for c in cands:
        if c in df.columns:
            return c
    raise ValueError(f"None of columns found: {cands}")


def parse_attrs(text):
    result={}
    for item in str(text).split(";"):
        if "=" in item:
            k,v=item.split("=",1); result[k]=v
    return result


def norm_contig(value,prefix):
    x=str(value)
    if "#" in x:
        x=x.split("#")[-2] if len(x.split("#"))>=3 else x
    if prefix and x.startswith(prefix):
        x=x[len(prefix):].lstrip("_")
    return x.lower()


def main():
    a=parse_args(); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    feat=pd.read_csv(a.features,sep="\t",dtype=str)
    if "community" not in feat.columns:
        raise ValueError("Feature metadata requires community column.")
    chrom_col=pick_col(feat,["chrom","CHROM","contig","reference_contig"])
    pos_col=pick_col(feat,["pos","POS","position","reference_position"])
    feat[pos_col]=pd.to_numeric(feat[pos_col],errors="coerce")
    feat=feat.dropna(subset=[pos_col]).copy()
    feat[pos_col]=feat[pos_col].astype(int)
    feat["_norm_chrom"]=feat[chrom_col].map(lambda x:norm_contig(x,a.reference_prefix))

    gff=pd.read_csv(a.gff3,sep="\t",comment="#",header=None,
                    names=["seqid","source","type","start","end","score","strand","phase","attributes"],
                    dtype={"seqid":str,"type":str,"attributes":str})
    gff=gff[gff["type"].isin(["gene","CDS","exon","intron","five_prime_UTR","three_prime_UTR"])].copy()
    gff["start"]=pd.to_numeric(gff["start"],errors="coerce")
    gff["end"]=pd.to_numeric(gff["end"],errors="coerce")
    gff=gff.dropna(subset=["start","end"])
    gff["_norm_chrom"]=gff["seqid"].map(lambda x:norm_contig(x,a.reference_prefix))
    gff["_attrs"]=gff["attributes"].map(parse_attrs)
    gff["gene_id"]=gff["_attrs"].map(lambda d:d.get("ID") or d.get("Parent") or d.get("locus_tag") or "")
    gff["product"]=gff["_attrs"].map(lambda d:d.get("product",""))

    priority={"CDS":1,"five_prime_UTR":2,"three_prime_UTR":2,"exon":3,"intron":4,"gene":5}
    hits=[]
    for chrom,sub in feat.groupby("_norm_chrom"):
        gf=gff[gff["_norm_chrom"]==chrom]
        if gf.empty:
            for idx,row in sub.iterrows():
                hits.append({**row.to_dict(),"region":"intergenic","gene_id":"","product":""})
            continue
        for idx,row in sub.iterrows():
            pos=row[pos_col]
            ov=gf[(gf["start"]<=pos)&(gf["end"]>=pos)].copy()
            if ov.empty:
                region="intergenic"; gene_id=""; product=""
            else:
                ov["_priority"]=ov["type"].map(priority).fillna(99)
                best=ov.sort_values("_priority").iloc[0]
                region=best["type"]; gene_id=best["gene_id"]; product=best["product"]
            hits.append({**row.to_dict(),"region":region,"gene_id":gene_id,"product":product})

    hit_df=pd.DataFrame(hits)
    drop_cols=[c for c in ["_norm_chrom"] if c in hit_df.columns]
    hit_df.drop(columns=drop_cols,inplace=True,errors="ignore")
    hit_df.to_csv(out/"community_snp_gene_hits.tsv",sep="\t",index=False)

    region=(hit_df.groupby(["community","region"]).size().rename("snp_count").reset_index())
    region["community_total"]=region.groupby("community")["snp_count"].transform("sum")
    region["percentage"]=region["snp_count"]/region["community_total"]*100
    region.to_csv(out/"community_region_summary.tsv",sep="\t",index=False)

    gene_hits=hit_df[hit_df["gene_id"].astype(str)!=""].copy()
    genes=(gene_hits.groupby(["community","gene_id","product"]).size().rename("snp_count").reset_index())
    genes.to_csv(out/"community_gene_summary.tsv",sep="\t",index=False)
    print(f"Annotated {len(hit_df)} SNP features.")


if __name__=="__main__":
    try: main()
    except Exception as e:
        print(f"[ERROR] {e}",file=sys.stderr); raise
