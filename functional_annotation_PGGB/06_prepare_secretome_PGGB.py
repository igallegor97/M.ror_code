#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def fasta(path):
    d={}; h=None; s=[]
    for raw in open(path):
        x=raw.strip()
        if not x: continue
        if x.startswith(">"):
            if h is not None: d[h]="".join(s)
            h=x[1:].split()[0]; s=[]
        else: s.append(x)
    if h is not None: d[h]="".join(s)
    return d

p=argparse.ArgumentParser()
p.add_argument("--proteins",required=True)
p.add_argument("--signalp-results",required=True)
p.add_argument("--deeptmhmm",required=True)
p.add_argument("--output-dir",required=True)
a=p.parse_args()
out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
seqs=fasta(a.proteins)
sp=pd.read_csv(a.signalp_results,sep="\t",comment="#")
sp.columns=[c.strip() for c in sp.columns]
idc="ID" if "ID" in sp.columns else sp.columns[0]
pc="Prediction" if "Prediction" in sp.columns else sp.columns[1]
sp=sp[[idc,pc]].rename(columns={idc:"protein_id",pc:"signalp_prediction"})
sp["signal_peptide"]=sp.signalp_prediction.astype(str).eq("SP")
tm=pd.read_csv(a.deeptmhmm,sep="\t")
x=sp.merge(tm,on="protein_id",how="left",validate="one_to_one")
x["transmembrane_domains"]=pd.to_numeric(x.transmembrane_domains,errors="coerce").fillna(0).astype(int)
x["secreted_candidate"]=x.signal_peptide & (x.transmembrane_domains==0)
x.to_csv(out/"secretome_classification.tsv",sep="\t",index=False)
ids=set(x.loc[x.secreted_candidate,"protein_id"])
with open(out/"MrorC26_secreted_candidates.fa","w") as h:
    for pid in sorted(ids):
        if pid not in seqs: continue
        h.write(f">{pid}\n")
        for i in range(0,len(seqs[pid]),80): h.write(seqs[pid][i:i+80]+"\n")
print(x.secreted_candidate.value_counts(dropna=False).to_string())
