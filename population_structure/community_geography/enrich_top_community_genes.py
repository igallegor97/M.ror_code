#!/usr/bin/env python3
"""
enrich_top_community_genes.py

Performs exploratory Fisher-exact enrichment for GO, PFAM and CAZyme terms
among genes hit by SNPs in top geographic communities. Also summarizes
secreted proteins and effector candidates.

Requires a curated gene annotation table with gene_id.
"""

from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact


def bh(p):
    s=pd.Series(p,dtype=float); order=s.sort_values().index; n=len(s)
    ranked=s.loc[order].values; adj=ranked*n/(range(1,n+1))
    import numpy as np
    adj=np.minimum.accumulate(adj[::-1])[::-1]; adj=np.minimum(adj,1)
    out=pd.Series(index=order,data=adj); return out.reindex(s.index)


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument("--gene-summary",required=True)
    p.add_argument("--annotations",required=True)
    p.add_argument("--top-communities",required=True)
    p.add_argument("--output-dir",required=True)
    return p.parse_args()


def split_terms(v):
    if pd.isna(v) or str(v).strip()=="":
        return []
    return [x.strip() for x in str(v).replace(",",";").split(";") if x.strip()]


def main():
    a=parse_args(); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    genes=pd.read_csv(a.gene_summary,sep="\t",dtype=str)
    ann=pd.read_csv(a.annotations,sep="\t",dtype=str).fillna("")
    top=pd.read_csv(a.top_communities,sep="\t",dtype=str)
    top_set=set(top["community"])
    hit=set(genes.loc[genes["community"].isin(top_set),"gene_id"])
    background=set(genes["gene_id"])
    ann=ann[ann["gene_id"].isin(background)].copy()

    enrichment=[]
    for category in ["GO","PFAM","CAZyme"]:
        if category not in ann.columns: continue
        term_to_genes={}
        for _,r in ann.iterrows():
            for term in split_terms(r[category]):
                term_to_genes.setdefault(term,set()).add(r["gene_id"])
        for term,term_genes in term_to_genes.items():
            a1=len(hit & term_genes)
            b1=len(hit - term_genes)
            c1=len((background-hit) & term_genes)
            d1=len((background-hit)-term_genes)
            odds,pv=fisher_exact([[a1,b1],[c1,d1]],alternative="greater")
            enrichment.append({"category":category,"term":term,"top_genes_with_term":a1,
                               "background_genes_with_term":len(term_genes),
                               "odds_ratio":odds,"p_value":pv})
    enr=pd.DataFrame(enrichment)
    if not enr.empty:
        enr["q_value"]=enr.groupby("category")["p_value"].transform(lambda x:bh(x).values)
        enr.sort_values(["q_value","p_value"]).to_csv(out/"functional_enrichment.tsv",sep="\t",index=False)
    else:
        pd.DataFrame(columns=["category","term","p_value","q_value"]).to_csv(out/"functional_enrichment.tsv",sep="\t",index=False)

    selected=ann[ann["gene_id"].isin(hit)].copy()
    selected.to_csv(out/"top_community_annotated_genes.tsv",sep="\t",index=False)
    summary=[]
    for col in ["secreted","effector_candidate"]:
        if col in selected.columns:
            truth=selected[col].astype(str).str.upper().isin(["TRUE","YES","1"])
            summary.append({"metric":col,"count":int(truth.sum()),"total_annotated_genes":len(selected)})
    pd.DataFrame(summary).to_csv(out/"top_community_candidate_summary.tsv",sep="\t",index=False)
    print(f"Top communities: {','.join(sorted(top_set))}; genes: {len(hit)}")


if __name__=="__main__":
    try: main()
    except Exception as e:
        print(f"[ERROR] {e}",file=sys.stderr); raise
