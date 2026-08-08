#!/usr/bin/env python3
import argparse, csv, re, sys
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument("input"); p.add_argument("output"); p.add_argument("id_map"); p.add_argument("stats")
a=p.parse_args()
seen=set(); n=aa=invalid=0; lengths=[]
allowed=set("ABCDEFGHIKLMNPQRSTVWXYZ*UOJBZ")
Path(a.output).parent.mkdir(parents=True,exist_ok=True)
Path(a.id_map).parent.mkdir(parents=True,exist_ok=True)
def records(path):
    head=None; seq=[]
    with open(path,encoding="utf-8") as fh:
        for raw in fh:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                if head is not None: yield head,"".join(seq)
                head=line[1:]; seq=[]
            else:
                if head is None: raise ValueError("Secuencia antes del primer encabezado")
                seq.append(re.sub(r"\s+","",line).upper())
    if head is not None: yield head,"".join(seq)
try:
  with open(a.output,"w",newline="\n") as out, open(a.id_map,"w",newline="",encoding="utf-8") as mp:
    w=csv.writer(mp,delimiter="\t"); w.writerow(["protein_id","original_header"])
    for head,seq in records(a.input):
      pid=head.split()[0]
      if not pid or re.search(r"[^A-Za-z0-9_.:|+-]",pid): raise ValueError(f"ID no admitido: {pid!r}")
      if pid in seen: raise ValueError(f"ID duplicado: {pid}")
      if not seq: raise ValueError(f"Secuencia vacía: {pid}")
      seen.add(pid); n+=1; aa+=len(seq); lengths.append(len(seq)); invalid+=sum(c not in allowed for c in seq)
      out.write(f">{pid}\n"); out.write("\n".join(seq[i:i+60] for i in range(0,len(seq),60))+"\n"); w.writerow([pid,head])
  if not n: raise ValueError("FASTA sin secuencias")
  with open(a.stats,"w",newline="",encoding="utf-8") as st:
    w=csv.writer(st,delimiter="\t"); w.writerow(["proteins","total_aa","min_length","mean_length","max_length","invalid_residues"])
    w.writerow([n,aa,min(lengths),f"{aa/n:.2f}",max(lengths),invalid])
  if invalid: raise ValueError(f"Se encontraron {invalid} residuos no admitidos")
except Exception as e:
  print(f"ERROR: {e}",file=sys.stderr); sys.exit(1)

