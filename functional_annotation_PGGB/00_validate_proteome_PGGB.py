#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, sys
from collections import Counter
from pathlib import Path
import pandas as pd

def read_fasta(path):
    h=None; s=[]
    for n,raw in enumerate(open(path),1):
        line=raw.strip()
        if not line: continue
        if line.startswith(">"):
            if h is not None: yield h,"".join(s)
            h=line[1:].split()[0]; s=[]
        else:
            if h is None: raise ValueError(f"Sequence before header at line {n}")
            s.append(line)
    if h is not None: yield h,"".join(s)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fasta",required=True)
    p.add_argument("--output-dir",required=True)
    p.add_argument("--min-length",type=int,default=10)
    a=p.parse_args()
    out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    records=list(read_fasta(a.fasta))
    ids=[r[0] for r in records]
    dup=sorted(k for k,v in Counter(ids).items() if v>1)
    allowed=set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    rows=[]; cleaned=[]
    for pid,seq in records:
        seq=seq.upper().rstrip("*")
        invalid="".join(sorted(set(seq)-allowed))
        rows.append({
            "protein_id":pid,"length_aa":len(seq),
            "internal_stop":"*" in seq,
            "invalid_residues":invalid,
            "sha256":hashlib.sha256(seq.encode()).hexdigest(),
            "passes_min_length":len(seq)>=a.min_length
        })
        if len(seq)>=a.min_length and not invalid and "*" not in seq:
            cleaned.append((pid,seq))
    df=pd.DataFrame(rows)
    df.to_csv(out/"proteome_qc_per_protein.tsv",sep="\t",index=False)
    summary=pd.DataFrame([
        ["proteins_total",len(records)],["unique_ids",len(set(ids))],
        ["duplicate_ids",len(dup)],
        ["short_proteins",int((~df.passes_min_length).sum())],
        ["internal_stops",int(df.internal_stop.sum())],
        ["invalid_residues",int(df.invalid_residues.ne("").sum())],
        ["cleaned_proteins",len(cleaned)]
    ],columns=["metric","value"])
    summary.to_csv(out/"proteome_qc_summary.tsv",sep="\t",index=False)
    (out/"duplicate_ids.txt").write_text("\n".join(dup)+("\n" if dup else ""))
    with open(out/"MrorC26_proteins.cleaned.fa","w") as handle:
        for pid,seq in cleaned:
            handle.write(f">{pid}\n")
            for i in range(0,len(seq),80): handle.write(seq[i:i+80]+"\n")
    print(summary.to_string(index=False))
    if dup: raise RuntimeError("Duplicate FASTA IDs detected.")

if __name__=="__main__":
    try: main()
    except Exception as e:
        print(f"[ERROR] {e}",file=sys.stderr); raise
