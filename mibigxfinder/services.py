# src/mibig_finder/services.py
"""Pure functions. CLI and future FastAPI both call these."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from mibigxfinder import chem, spectra
from mibigxfinder.config import BgcSearchConfig, SpectraSearchConfig, StructureSearchConfig


# ---------------------------------------------------------------------------
# Spectra search
# ---------------------------------------------------------------------------

def search_spectra(cfg: SpectraSearchConfig) -> pd.DataFrame:
    """Load query + reference MGFs, score, return top-N matches above threshold."""
    if cfg.clean:
        cleaned_path = cfg.out_dir / "cleaned_queries.mgf"
        cfg.out_dir.mkdir(parents=True, exist_ok=True)
        queries = spectra.cleaning_spectra(cfg.query, str(cleaned_path))
    else:
        queries = spectra.load_spectra(cfg.query)
    references = spectra.load_spectra(cfg.database) 

    scores = spectra.score(
        queries=queries,
        references=references,
        similarity=cfg.similarity,
        tolerance=cfg.tolerance,
    )

    score_arr = scores.to_array("score")      # (n_references, n_queries)
    match_arr = scores.to_array("matches")

    rows = []
    for i, q in enumerate(queries):
        q_id = q.metadata.get("scans") or q.metadata.get("feature_id") or i

        col = score_arr[:, i]
        order = np.argsort(col)[::-1]

        kept = 0
        for ref_idx in order:
            score_val = float(col[ref_idx])
            if score_val < cfg.sim_threshold:
                break
            kept += 1
            ref = references[ref_idx]
            rows.append({
                "feature_id": q_id,
                "library_id": ref.metadata.get("spectrum_id") or ref.metadata.get("spectrum_idx"),
                "inchikey": ref.metadata.get("inchikey"),
                "spectral similarity": score_val,
                "matched_peaks": int(match_arr[ref_idx, i]),
                "rank": kept,
            })
            if kept >= cfg.top_n:
                break

    df = pd.DataFrame(rows)

    sep = "\t" if cfg.bgc_links.suffix.lower() == ".tsv" else ","
    links = pd.read_csv(cfg.bgc_links, sep=sep)[["spectrum_idx", "bgc_id", "bgc_compound"]]
    links["spectrum_idx"] = links["spectrum_idx"].astype(str)
    df["library_id"] = df["library_id"].astype(str)

    df = df.merge(
        links,
        left_on="library_id", right_on="spectrum_idx",
        how="left",
    ).drop(columns=["spectrum_idx"])

    df = df.sort_values(["feature_id", "spectral similarity"], ascending=[True, False])
    df["rank"] = df.groupby("feature_id").cumcount() + 1
    df = df[["feature_id", "library_id", "bgc_id", "bgc_compound",
            "inchikey", "spectral similarity", "matched_peaks", "rank"]]    
    _maybe_write(df, cfg.out_dir, name="spectra_matches")
    return df


# ---------------------------------------------------------------------------
# Structure search
# ---------------------------------------------------------------------------

def search_structures(cfg: StructureSearchConfig) -> pd.DataFrame:
    """Load query SMILES, score against reference library (SMILES pulled from MGF metadata)."""
    queries = _read_smiles(cfg.query)                       # [(id, smiles), ...]
    refs = _refs_from_mgf(cfg.database)                     # [(library_id, bgc_id, smiles), ...]

    # precompute reference fingerprints per fp kind
    ref_fps = {
        kind: [(lib_id, bgc, chem.compute_fp(smi, kind)) for lib_id, bgc, smi in refs]
        for kind in cfg.fps
    }

    rows = []
    for qid, qsmi in queries:
        if qsmi.lower() == "smiles":
            continue
        qsmi_c = chem.canonical_smiles(qsmi)
        if qsmi_c is None:
            continue
        qkey = chem.inchikey14(qsmi_c)

        per_kind = {kind: chem.compute_fp(qsmi_c, kind) for kind in cfg.fps}

        for i, (lib_id, bgc, rsmi) in enumerate(refs):
            sims = {}
            for kind in cfg.fps:
                q_fp = per_kind[kind]
                r_fp = ref_fps[kind][i][2]
                if q_fp is None or r_fp is None:
                    continue
                sims[kind] = chem.tanimoto(q_fp, r_fp)
            if not sims:
                continue
            best = max(sims.values())
            if best < cfg.similarity_threshold:
                continue
            rows.append({
                "query_id": qid,
                "query_inchikey14": qkey,
                "library_id": lib_id,
                "bgc_match": bgc,
                "best_tanimoto": best,
                **{f"tanimoto_{k}": v for k, v in sims.items()},
            })

        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values(["query_id", "best_tanimoto"], ascending=[True, False])
        _maybe_write(df, cfg.out_dir, name="structure_matches")
    return df


# ---------------------------------------------------------------------------
# BGC search
# ---------------------------------------------------------------------------

def search_bgcs(cfg: BgcSearchConfig) -> pd.DataFrame:
    """Rank MIBiG BGCs by BiG-SCAPE distance to each query region; keep the closest ones with spectra.

    One row per (query, MIBiG BGC, linked spectrum). ``mibig_rank`` is the rank among
    all MIBiG BGCs, so a gap to 1 means closer BGCs without spectra exist.
    Writes ``bgc_matches.{csv,parquet}`` and the per-query ``bgc_summary.{csv,parquet}``.
    """
    from mibigxfinder import bgc            # needs BiG-SCAPE, so only imported here

    manifest = cfg.reference.with_suffix(".json")
    if manifest.exists() and bgc.sha256(cfg.pfam) != json.loads(manifest.read_text())["pfam_sha256"]:
        raise ValueError(
            f"{cfg.pfam} is not the Pfam file {cfg.reference.name} was built with (see {manifest.name}); "
            "queries and reference must use the same Pfam release"
        )
    references = pd.read_parquet(cfg.reference)
    sep = "\t" if cfg.bgc_links.suffix.lower() == ".tsv" else ","
    links = pd.read_csv(cfg.bgc_links, sep=sep)
    missing = sorted(set(links["bgc_id"]) - set(references["bgc_id"]))
    if missing:
        warnings.warn(f"{len(missing)} linked BGCs are not in {cfg.reference.name} and can't be found, e.g. {missing[:3]}")

    hits = bgc.score(bgc.annotate(cfg.query, cfg.pfam, cores=cfg.cores), references)
    hits = hits.rename(columns={"reference_id": "bgc_id"})
    hits["mibig_rank"] = hits.groupby("query_id").cumcount() + 1   # hits are sorted by distance
    hits = hits[(hits["distance"] < cfg.max_distance) & hits["bgc_id"].isin(links["bgc_id"])]
    hits = hits.groupby("query_id").head(cfg.top_n)

    df = hits.merge(links, on="bgc_id")
    _maybe_write(df, cfg.out_dir, name="bgc_matches")
    _maybe_write(summarize_bgc_matches(df, [f.stem for f in cfg.query]), cfg.out_dir, name="bgc_summary")
    return df


def summarize_bgc_matches(matches: pd.DataFrame, query_ids: list[str]) -> pd.DataFrame:
    """One row per query: its closest MIBiG BGC with spectra and how many BGCs/spectra it got.

    Queries without matches get a row too, with ``n_bgcs`` 0.
    """
    rows = []
    for query_id in query_ids:
        m = matches[matches["query_id"] == query_id].sort_values(["distance", "mibig_rank"], kind="stable")
        row = {"query_id": query_id, "closest_bgc": None, "distance": float("nan"),
               "mibig_rank": pd.NA, "compound": None, "n_bgcs": 0, "n_spectra": 0}
        if len(m):
            best = m.iloc[0]
            row.update(closest_bgc=best["bgc_id"], distance=best["distance"], mibig_rank=best["mibig_rank"],
                       n_bgcs=m["bgc_id"].nunique(),
                       n_spectra=m["spectrum_idx"].nunique() if "spectrum_idx" in m else len(m))
            if "bgc_compound" in m:
                names = list(dict.fromkeys(m.loc[m["bgc_id"] == best["bgc_id"], "bgc_compound"].dropna()))
                if names:
                    row["compound"] = names[0] + (f" (+{len(names) - 1})" if len(names) > 1 else "")
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary["mibig_rank"] = summary["mibig_rank"].astype("Int64")
    if summary["compound"].isna().all():
        summary = summary.drop(columns="compound")      # links file without compound names
    return summary


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# src/mibig_finder/services.py

import pandas as pd
def _read_smiles(path: Path) -> list[tuple[str, str]]:
    with path.open() as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    if lines and lines[0].lower().startswith("smiles"):
        lines = lines[1:]
    if lines and ("," in lines[0] or "\t" in lines[0]):
        sep = "," if "," in lines[0] else "\t"
        rows = [l.split(sep) for l in lines]
        if len(rows[0]) >= 2:
            return [(r[0].strip(), r[1].strip()) for r in rows]
    return [(smi, smi) for smi in lines]

def _refs_from_mgf(path: Path) -> list[tuple[str, str | None, str]]:
    """Pull (library_id, bgc_id, smiles) from each reference spectrum's metadata."""
    out = []
    for s in spectra.load_spectra(path):
        smi = s.metadata.get("smiles") 
        if not smi:
            continue
        lib_id = s.metadata.get("spectrum_id") or s.metadata.get("title") or ""
        bgc = s.metadata.get("bgc_match")
        out.append((lib_id, bgc, smi))
    return out

def _maybe_write(df: pd.DataFrame, out_dir: Path | None, name: str) -> None:
    if out_dir is None:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"{name}.csv", index=False)
    df.to_parquet(out_dir / f"{name}.parquet", index=False)