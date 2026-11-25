# mibigxfinder/bgc.py
"""BGC distances with BiG-SCAPE 2's metric, without its database.

``annotate`` turns region .gbk files into a domain table (one row per CDS and
Pfam domain) using BiG-SCAPE's own parsing and HMMER code. ``score`` rebuilds
BiG-SCAPE records from such tables and runs BiG-SCAPE's own pair scoring.
Settings are the ``bigscape query`` defaults: region records, glocal
alignment, legacy extend strategy, mix weights.

Needs BiG-SCAPE 2 (bioconda) and a Pfam-A.hmm; tables are only comparable if
built with the same Pfam release and BiG-SCAPE version.
"""

from __future__ import annotations

import hashlib
import json
import warnings
from datetime import date
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

import pandas as pd

# BiG-SCAPE calls warnings.filterwarnings("ignore") when its profiler module is imported,
# which would silence warnings in the whole process; catch_warnings undoes that afterwards.
with warnings.catch_warnings():
    import big_scape
    from big_scape.cli.config import BigscapeConfig
    from big_scape.comparison.workflow import calculate_scores_pair
    from big_scape.enums import ALIGNMENT_MODE, EXTEND_STRATEGY, SOURCE_TYPE
    from big_scape.genbank import CDS, GBK, Region
    from big_scape.hmm import HMMer, HSP, HSPAlignment

CONFIG = Path(big_scape.__file__).parent / "config.yml"
ALIGNMENT, EXTEND, WEIGHTS = ALIGNMENT_MODE.GLOCAL, EXTEND_STRATEGY.LEGACY, "mix"

COLUMNS = [
    "bgc_id", "region_start", "region_stop", "contig_edge",
    "orf_num", "cds_start", "cds_stop", "strand", "gene_kind",
    "domain", "score", "env_start", "env_stop", "alignment",
]

# GBK.parse only reads these run settings
_PARSE_RUN = {
    "classify": None, "legacy_weights": False, "force_gbk": False,
    "include_categories": None, "exclude_categories": None,
}


def annotate(
    gbk_files: Iterable[str | Path], pfam_path: str | Path, *, mibig: bool = False, cores: int = 1
) -> pd.DataFrame:
    """Domain table for region .gbk files; ``bgc_id`` is the file name without suffix.

    Use ``mibig=True`` for MIBiG files: BiG-SCAPE then accepts regions without
    candidate cluster features, as it does for ``--mibig-version``.
    """
    BigscapeConfig.parse_config(CONFIG)
    source = SOURCE_TYPE.MIBIG if mibig else SOURCE_TYPE.QUERY
    gbks = [
        GBK.parse(Path(f), source, _PARSE_RUN, BigscapeConfig.CDS_OVERLAP_CUTOFF)
        for f in gbk_files
    ]
    genes = [cds for gbk in gbks for cds in gbk.genes]

    # same steps as bigscape's hmmsearch/hmmalign tasks (simple search on every platform)
    pfam_path = Path(pfam_path)
    if not HMMer.are_profiles_pressed(pfam_path):
        HMMer.press(pfam_path)
    HMMer.init(pfam_path)
    HMMer.hmmsearch_simple(genes, BigscapeConfig.DOMAIN_OVERLAP_CUTOFF, cores)
    HMMer.unload()
    for cds in genes:
        cds.hsps = sorted(cds.hsps)
    HMMer.init(pfam_path, False)
    HMMer.align_simple([hsp for cds in genes for hsp in cds.hsps])
    HMMer.unload()

    rows = []
    for gbk in gbks:
        r = gbk.region
        region = (gbk.path.stem, r.nt_start, r.nt_stop, bool(r.contig_edge))
        for cds in gbk.genes:
            gene = (cds.orf_num, cds.nt_start, cds.nt_stop, cds.strand, cds.gene_kind)
            for hsp in cds.hsps or [None]:
                domain = (
                    (hsp.domain, hsp.score, hsp.env_start, hsp.env_stop, hsp.alignment.align_string)
                    if hsp is not None else (None,) * 5
                )
                rows.append(region + gene + domain)
    return pd.DataFrame(rows, columns=COLUMNS)


def _records(table: pd.DataFrame) -> dict[str, Region]:
    """Rebuild minimal BiG-SCAPE region records, like bigscape's fetch_records_from_database."""
    records: dict[str, Region] = {}
    cds = None
    for row in table[COLUMNS].itertuples(index=False):
        if row.bgc_id not in records:
            gbk = GBK(Path(row.bgc_id), "", SOURCE_TYPE.REFERENCE)
            records[row.bgc_id] = Region(
                gbk, 0, int(row.region_start), int(row.region_stop), bool(row.contig_edge), ""
            )
            cds = None
        gbk = records[row.bgc_id].parent_gbk
        if cds is None or cds.orf_num != row.orf_num:
            cds = CDS(int(row.cds_start), int(row.cds_stop))
            cds.orf_num, cds.strand = int(row.orf_num), int(row.strand)
            cds.gene_kind = row.gene_kind if isinstance(row.gene_kind, str) else None
            gbk.genes.append(cds)
        if isinstance(row.domain, str):
            hsp = HSP(cds, row.domain, float(row.score), int(row.env_start), int(row.env_stop))
            hsp.alignment = HSPAlignment(hsp, row.alignment)
            cds.hsps.append(hsp)
    return records


def score(queries: pd.DataFrame, references: pd.DataFrame) -> pd.DataFrame:
    """BiG-SCAPE distances of every query BGC to every reference BGC, closest first."""
    cols = ["query_id", "reference_id", "distance", "jaccard", "adjacency", "dss"]
    query_recs, ref_recs = _records(queries), _records(references)
    names = {}
    for i, (name, rec) in enumerate([*query_recs.items(), *ref_recs.items()]):
        rec._db_id = i          # scoring indexes records by this id
        names[i] = name
    pairs = [(q, r) for q in query_recs.values() for r in ref_recs.values()]
    if not pairs:
        return pd.DataFrame(columns=cols)

    results = calculate_scores_pair((pairs, ALIGNMENT, EXTEND, 0, WEIGHTS, CONFIG))
    df = pd.DataFrame(
        [(names[a], names[b], dist, jc, ai, dss) for a, b, dist, jc, ai, dss, *_ in results],
        columns=cols,
    )
    return df.sort_values(["query_id", "distance", "reference_id"], ignore_index=True)


def build_reference(
    out_dir: str | Path,
    pfam_path: str | Path,
    *,
    mibig_version: str = "5.0",
    gbk_dir: str | Path | None = None,
    cores: int = 1,
    batch: int = 100,
) -> Path:
    """Annotate MIBiG .gbk files into ``<out_dir>/mibig_<version>.parquet`` plus a JSON manifest.

    Without ``gbk_dir`` the antiSMASH-processed MIBiG release is downloaded into
    ``<out_dir>/MIBiG/`` by BiG-SCAPE's own downloader (reused if already there).
    """
    out_dir = Path(out_dir)
    if gbk_dir is None:
        from big_scape.file_input.load_files import get_mibig
        print(f"getting antiSMASH-processed MIBiG {mibig_version} ...")
        gbk_dir = get_mibig(mibig_version, out_dir)
    files = sorted(Path(gbk_dir).glob("*.gbk"))
    if not files:
        raise FileNotFoundError(f"no .gbk files in {gbk_dir}")

    parts = []
    for i in range(0, len(files), batch):   # batches only to report progress
        parts.append(annotate(files[i:i + batch], pfam_path, mibig=True, cores=cores))
        print(f"annotated {min(i + batch, len(files))}/{len(files)} BGCs")
    table = pd.concat(parts, ignore_index=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"mibig_{mibig_version}.parquet"
    table.to_parquet(out, index=False)
    manifest = {
        "mibig_version": mibig_version,
        "gbk_dir": Path(gbk_dir).name,
        "n_bgcs": int(table["bgc_id"].nunique()),
        "n_domains": int(table["domain"].notna().sum()),
        "pfam_sha256": sha256(pfam_path),
        "bigscape_version": version("big-scape"),
        "settings": {
            "record_type": "region", "alignment_mode": ALIGNMENT.value,
            "extend_strategy": EXTEND.value, "weights": WEIGHTS,
        },
        "created": date.today().isoformat(),
    }
    out.with_suffix(".json").write_text(json.dumps(manifest, indent=2))
    return out


def sha256(path: str | Path) -> str:
    """Checksum of a file, e.g. to check that queries use the reference's Pfam release."""
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()
