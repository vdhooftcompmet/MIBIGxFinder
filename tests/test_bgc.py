"""BGC search. Needs BiG-SCAPE; the end-to-end test also needs data/pfam/Pfam-A.hmm.

Fixtures: NC_003888.3.region002 (antiSMASH region of S. coelicolor A3(2)) as query, its
domain table, 23 MIBiG 4.0 BGCs as reference, and the distances BiG-SCAPE 2.0.3 query
mode gave for them (`bgc_distances_bigscape.tsv`).
"""

from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("big_scape")

from mibigxfinder import bgc  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"
PFAM = Path(__file__).parents[1] / "data/pfam/Pfam-A.hmm"
QUERY = FIXTURES / "NC_003888.3.region002.gbk"
REFERENCE = FIXTURES / "bgc_reference_small.parquet"


def test_score_reproduces_bigscape_query_mode():
    query = pd.read_parquet(FIXTURES / "NC_003888.3.region002.parquet")
    hits = bgc.score(query, pd.read_parquet(REFERENCE))
    expected = pd.read_csv(FIXTURES / "bgc_distances_bigscape.tsv", sep="\t")

    merged = expected.merge(hits, left_on="bgc_id", right_on="reference_id", validate="1:1")
    assert len(merged) == len(expected) == len(hits)
    assert (merged["distance_x"] - merged["distance_y"]).abs().max() < 1e-12


@pytest.mark.skipif(not PFAM.exists(), reason="needs data/pfam/Pfam-A.hmm")
def test_bgc_search_returns_linked_spectra(tmp_path):
    from mibigxfinder.config import BgcSearchConfig
    from mibigxfinder.services import search_bgcs

    cfg = BgcSearchConfig(
        query=[QUERY], reference=REFERENCE, bgc_links=FIXTURES / "bgc_links.csv",
        pfam=PFAM, top_n=2, out_dir=tmp_path,
    )
    with pytest.warns(UserWarning, match="linked BGCs are not in"):   # BGC0000005, BGC0003163, ...
        df = search_bgcs(cfg)

    assert df["bgc_id"].tolist() == ["BGC0000664", "BGC0000649"]
    assert df["spectrum_idx"].tolist() == ["CCMSLIB00000079400", "CCMSLIB00000079401"]
    assert df["mibig_rank"].tolist() == [1, 2]
    assert (tmp_path / "bgc_matches.csv").exists()
