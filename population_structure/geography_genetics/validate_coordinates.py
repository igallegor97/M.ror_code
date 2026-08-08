#!/usr/bin/env python3
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd

ALLOWED_PRECISION={"locality","municipality","state","country","missing"}

def args():
    p=argparse.ArgumentParser()
    p.add_argument("--coordinates",required=True)
    p.add_argument("--metadata",required=True)
    p.add_argument("--output-dir",required=True)
    return p.parse_args()

def main():
    a=args()
    out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    coords=pd.read_csv(a.coordinates,sep="\t",dtype=str,keep_default_na=False)
    meta=pd.read_csv(a.metadata,sep="\t",dtype=str)
    required={"sample_id","country","region_state","latitude","longitude",
              "coordinate_precision","coordinate_source","notes"}
    missing=required-set(coords.columns)
    if missing: raise ValueError(f"Missing columns: {sorted(missing)}")
    if coords.sample_id.duplicated().any():
        raise ValueError("Duplicated sample IDs.")
    expected=set(meta.sample_id); observed=set(coords.sample_id)
    if expected!=observed:
        raise ValueError(f"Sample mismatch. Missing={sorted(expected-observed)}; extra={sorted(observed-expected)}")
    bad=set(coords.coordinate_precision)-ALLOWED_PRECISION
    if bad: raise ValueError(f"Invalid precision values: {sorted(bad)}")

    issues=[]; valid=[]
    for _,r in coords.iterrows():
        sid=r.sample_id; precision=r.coordinate_precision
        lat=r.latitude.strip(); lon=r.longitude.strip()
        if precision=="missing":
            if lat or lon: issues.append([sid,"missing_precision_but_coordinates_present"])
            continue
        if not lat or not lon:
            issues.append([sid,"coordinates_missing_for_nonmissing_precision"]); continue
        try:
            latf=float(lat); lonf=float(lon)
        except ValueError:
            issues.append([sid,"coordinates_not_numeric"]); continue
        if not -90<=latf<=90:
            issues.append([sid,"latitude_out_of_range"]); continue
        if not -180<=lonf<=180:
            issues.append([sid,"longitude_out_of_range"]); continue
        if not r.coordinate_source.strip():
            issues.append([sid,"coordinate_source_missing"])
        valid.append(r.to_dict())

    pd.DataFrame(issues,columns=["sample_id","issue"]).to_csv(
        out/"coordinate_validation_issues.tsv",sep="\t",index=False)
    pd.DataFrame(valid,columns=coords.columns).to_csv(
        out/"sample_coordinates_24G_validated.tsv",sep="\t",index=False)
    summary=pd.DataFrame([
        ["metadata_samples",len(expected)],
        ["coordinate_rows",len(coords)],
        ["valid_coordinate_rows",len(valid)],
        ["missing_coordinate_rows",(coords.coordinate_precision=="missing").sum()],
        ["validation_issues",len(issues)],
        ["locality_or_municipality_rows",coords.coordinate_precision.isin(["locality","municipality"]).sum()]
    ],columns=["metric","value"])
    summary.to_csv(out/"coordinate_validation_summary.tsv",sep="\t",index=False)
    print(summary.to_string(index=False))
    if len(valid)<3: raise RuntimeError("Fewer than three valid coordinates.")

if __name__=="__main__":
    try: main()
    except Exception as e:
        print(f"[ERROR] {e}",file=sys.stderr,flush=True); raise
