#!/usr/bin/env python3
import argparse,re,pandas as pd
p=argparse.ArgumentParser()
p.add_argument("--three-line",required=True)
p.add_argument("--output",required=True)
a=p.parse_args()
lines=[x.strip() for x in open(a.three_line) if x.strip()]
if len(lines)%3: raise ValueError("Expected DeepTMHMM 3-line format.")
rows=[]
for i in range(0,len(lines),3):
    pid=lines[i].lstrip(">").split()[0]
    seq,topo=lines[i+1],lines[i+2]
    if len(seq)!=len(topo): raise ValueError(f"Length mismatch: {pid}")
    runs=re.findall(r"[MB]+",topo)
    rows.append({"protein_id":pid,"transmembrane_domains":len(runs),
                 "has_transmembrane_domain":bool(runs),"topology":topo})
pd.DataFrame(rows).to_csv(a.output,sep="\t",index=False)
