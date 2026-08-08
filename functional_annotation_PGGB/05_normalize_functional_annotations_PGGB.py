#!/usr/bin/env python3
import argparse,csv,re
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument("--root",required=True); p.add_argument("--samples",required=True); a=p.parse_args(); root=Path(a.root)
samples=[]
with open(a.samples,encoding="utf-8") as f:
 for r in csv.DictReader(f,delimiter="\t"): samples.append(r["sample_id"])
norm=root/"05_normalized"; norm.mkdir(parents=True,exist_ok=True)
master=[]
for s in samples:
 fasta=root/"00_qc/fasta"/f"{s}.faa"; ids=[x[1:].split()[0] for x in fasta.read_text().splitlines() if x.startswith(">")]
 ann={i:{"sample_id":s,"protein_id":i,"eggnog_description":"","GO_terms":"","KEGG_ko":"","pfam_domains":"","cazy_families":"","signalp_prediction":""} for i in ids}
 eg=next(iter((root/"01_eggnog"/s).glob("*.emapper.annotations")),None)
 if eg:
  with eg.open(encoding="utf-8") as f:
   header=None
   for line in f:
    if line.startswith("##"): continue
    if line.startswith("#query"): header=line[1:].rstrip().split("\t"); continue
    if not header or line.startswith("#"): continue
    d=dict(zip(header,line.rstrip().split("\t"))); q=d.get("query")
    if q in ann:
     ann[q]["eggnog_description"]=d.get("Description",""); ann[q]["GO_terms"]=d.get("GOs",""); ann[q]["KEGG_ko"]=d.get("KEGG_ko","")
 pf=root/"02_pfam"/s/f"{s}.pfam.domtblout"
 if pf.exists():
  hits={}
  for line in pf.read_text().splitlines():
   if line.startswith("#") or not line.strip(): continue
   x=line.split(); hits.setdefault(x[3],set()).add(x[0])
  for q,v in hits.items():
   if q in ann: ann[q]["pfam_domains"]=";".join(sorted(v))
 dbdir=root/"03_dbcan"/s
 for candidate in [dbdir/"overview.tsv",dbdir/"overview.txt",dbdir/"CAZyme_annotation.tsv"]:
  if candidate.exists():
   with candidate.open(encoding="utf-8") as f:
    rows=list(csv.reader(f,delimiter="\t"))
   if rows:
    hdr=rows[0]
    for row in rows[1:]:
     if not row: continue
     q=row[0]; fam=set(re.findall(r"(?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?",";".join(row)))
     if q in ann: ann[q]["cazy_families"]=";".join(sorted(fam))
   break
 spdir=root/"04_signalp"/s
 for candidate in list(spdir.glob("*summary*"))+list(spdir.glob("*prediction*"))+list(spdir.glob("*.txt")):
  for line in candidate.read_text(errors="ignore").splitlines():
   if not line or line.startswith("#"): continue
   x=line.split();
   if x and x[0] in ann: ann[x[0]]["signalp_prediction"]=" ".join(x[1:])
 out=norm/f"{s}.functional.normalized.tsv"; fields=list(next(iter(ann.values())).keys())
 with out.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(ann.values())
 master.extend(ann.values())
out=root/"07_master/functional_annotation_master.tsv"; out.parent.mkdir(parents=True,exist_ok=True)
fields=["sample_id","protein_id","eggnog_description","GO_terms","KEGG_ko","pfam_domains","cazy_families","signalp_prediction"]
with out.open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(master)

