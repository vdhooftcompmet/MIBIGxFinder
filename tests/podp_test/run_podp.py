"""Runs MIBiG-Finder on MIBiG-Finder/tests/podp_test, to see how the Finder works.

Usage (from anywhere):
    python run_podp.py [out_dir]     # default out_dir: ./structure_csv

"""
import sys, warnings, logging
warnings.filterwarnings("ignore")
logging.getLogger("matchms").setLevel(logging.ERROR)
from pathlib import Path
import pandas as pd
import csv
import requests
import shutil
import zipfile

from mibig_finder.config import StructureSearchConfig
from mibig_finder.services import search_structures

import find_links

CONCEPT_RECORD_ID = "7497442"

def get_latest_zenodo_record() -> dict:
    """Resolve the concept DOI to the metadata of the newest version."""
    r = requests.get(f"https://zenodo.org/api/records/{CONCEPT_RECORD_ID}", timeout=30)
    r.raise_for_status()
    meta = r.json()
    latest_url = meta.get("links", {}).get("latest")
    if latest_url:
        r2 = requests.get(latest_url, timeout=30)
        r2.raise_for_status()
        return r2.json()
    return meta

HERE = Path(__file__).resolve().parent
fx = HERE 

RAW_DIR = HERE / "raw_json"
OUT_TSV = HERE / "mibig_smiles_links.tsv"

# check if the file with the structures exists
rows = []

zip_path = RAW_DIR / "database.zip"

if not RAW_DIR.exists() or not zip_path.exists():
    RAW_DIR.mkdir(exist_ok=True)    

    if zip_path.exists():
        print(f"Using existing {zip_path}")
    else:
        record = get_latest_zenodo_record()
        version = record.get("metadata", {}).get("version", "unknown")
        files = record.get("files", [])
        db_file = next((f for f in files if f["key"].endswith(".zip")), None)
        if db_file is None:
            raise RuntimeError("No .zip file in the Zenodo record -- check manually")
        url = db_file["links"]["self"]
        print(f"Downloading PoDP database ({version}) from {url}")
        with requests.get(url, stream=True, timeout=120) as resp:
            resp.raise_for_status()
            with open(zip_path, "wb") as f:
                shutil.copyfileobj(resp.raw, f)

    with zipfile.ZipFile(zip_path) as z:
            z.extractall(RAW_DIR)

for json_path in sorted(RAW_DIR.glob("*.json")):
    rows.extend(find_links.extract_rows(json_path))

with open(OUT_TSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=find_links.COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

n_projects = len({r["json_stem"] for r in rows})
print(f"Wrote {len(rows)} rows from {n_projects} projects to {OUT_TSV}")


out = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "podp_out"
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_colwidth", 40)

# Finder input: two columns (id, smiles) with a header. The id is the row number,
# so each hit can be joined back to its PoDP row (same SMILES occurs in several projects).
QUERY_TSV = HERE / "podp_structures.tsv"
with open(QUERY_TSV, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["id", "smiles"])
    for i, r in enumerate(rows):
        writer.writerow([f"r{i}", r["SMILES"]])
print(f"Wrote {len(rows)} SMILES to {QUERY_TSV}")

print("\n=== STRUCTURE SEARCH: on PODP data ===")
df = search_structures(StructureSearchConfig(query=QUERY_TSV, database=fx/"<LIBRARY>.mgf", out_dir=out/"structure_csv"))
# TODO: <LIBRARY>.mgf because I am not sure how to get the library yet

# compare with the real BGCs
# Join the Finder's bgc_match back to the PoDP MIBiG column through the query id (r0, r1, ...)
truth = pd.DataFrame(rows)
truth["query_id"] = [f"r{i}" for i in range(len(rows))] #  the ID in the query

if df.empty:
    hit_bgcs = pd.Series(dtype=object)
else:
    hit_bgcs = df.groupby("query_id")["bgc_match"].agg(lambda s: set(s.dropna()))

truth["hit_bgcs"] = truth["query_id"].map(hit_bgcs) # the set of BGCs the Finder returned for that query
truth["n_hits"] = truth["hit_bgcs"].map(lambda s: len(s) if isinstance(s, set) else 0) # how many different BGCs are in hit_bgcs
truth["podp_bgc_found"] = [ # True if the PoDP MIBiG value is in hit_bgcs
    isinstance(hits, set) and bgc in hits
    for hits, bgc in zip(truth["hit_bgcs"], truth["MIBiG"])
]

# How often does the PoDP BGC appear among the hits?
n_queries = len(truth)
n_with_hits = int((truth["n_hits"] > 0).sum())
n_found = int(truth["podp_bgc_found"].sum())
print("\n=== PoDP BGC among the Finder hits ===")
print(f"Queries:                 {n_queries}")
print(f"Queries with any hit:    {n_with_hits}")
print(f"PoDP BGC among the hits: {n_found} ({n_found / n_queries:.1%})")

COMPARISON_TSV = out / "podp_comparison.tsv"
out.mkdir(parents=True, exist_ok=True)
truth.assign(hit_bgcs=truth["hit_bgcs"].map(lambda s: ",".join(sorted(s)) if isinstance(s, set) else "")) \
     .to_csv(COMPARISON_TSV, sep="\t", index=False)
print(f"Wrote per-query comparison to {COMPARISON_TSV}") 