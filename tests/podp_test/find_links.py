import csv
import json
from pathlib import Path

DATA_DIR = Path("./podp_data")
RAW_DIR = DATA_DIR / "raw_json"
OUT_TSV = DATA_DIR / "mibig_smiles_links.tsv"

COLUMNS = [
    "json_stem",
    "GNPSMassIVE_ID",
    "molecular_network",
    "MIBiG",
    "match_type",
    "SMILES",
    "IUPAC",
    "MS2_scan",
    "link",
]


def extract_rows(json_path: Path) -> list[dict]:
    """Return one row per BGC_MS2_links entry that has a MIBiG number and SMILES."""
    with open(json_path) as f:
        data = json.load(f)

    project = data.get("metabolomics", {}).get("project", {}) or {}
    rows = []
    for entry in data.get("BGC_MS2_links", []) or []:
        smiles = entry.get("SMILES")
        if not smiles:
            continue

        bgc = entry.get("BGC_ID", {}) or {}
        if bgc.get("MIBiG_number"):
            mibig, match_type = bgc["MIBiG_number"], "exact"
        elif bgc.get("similar_MIBiG_number"):
            mibig, match_type = bgc["similar_MIBiG_number"], "similar"
        else:
            continue

        rows.append({
            "json_stem": json_path.stem,
            "GNPSMassIVE_ID": project.get("GNPSMassIVE_ID", ""),
            "molecular_network": project.get("molecular_network", ""),
            "MIBiG": mibig,
            "match_type": match_type,
            "SMILES": smiles,
            "IUPAC": entry.get("IUPAC", ""),
            "MS2_scan": entry.get("MS2_scan", ""),
            "link": entry.get("link", ""),
        })
    return rows