"""BGC search. Needs BiG-SCAPE.

Fixtures: NC_003888.3.region002 (antiSMASH region of S. coelicolor A3(2)) as query, its
domain table, 23 MIBiG 4.0 BGCs as reference, and the distances BiG-SCAPE 2.0.3 query
mode gave for them (`bgc_distances_bigscape.tsv`). All were made with `pfam_small.hmm`:
the 21 Pfam families BiG-SCAPE finds in the query region, copied unchanged from the Pfam
release downloaded on 2026-11-27 (the query gets the same domains as with the full Pfam).
"""

import shutil
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("big_scape")

from mibigxfinder import bgc  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"
QUERY = FIXTURES / "NC_003888.3.region002.gbk"
QUERY_TABLE = FIXTURES / "NC_003888.3.region002.parquet"
REFERENCE = FIXTURES / "bgc_reference_small.parquet"
BIGSCAPE = FIXTURES / "bgc_distances_bigscape.tsv"


@pytest.fixture(scope="module")
def pfam(tmp_path_factory):
    """The small Pfam in a temp folder, so its pressed .h3* files don't land in the repo."""
    path = tmp_path_factory.mktemp("pfam") / "pfam_small.hmm"
    shutil.copy(FIXTURES / "pfam_small.hmm", path)
    return path


def test_score_reproduces_bigscape_query_mode():
    hits = bgc.score(pd.read_parquet(QUERY_TABLE), pd.read_parquet(REFERENCE))
    expected = pd.read_csv(BIGSCAPE, sep="\t")

    merged = expected.merge(hits, left_on="bgc_id", right_on="reference_id", validate="1:1")
    assert len(merged) == len(expected) == len(hits)
    assert (merged["distance_x"] - merged["distance_y"]).abs().max() < 1e-12


def test_annotate_reproduces_query_table(pfam):
    pd.testing.assert_frame_equal(bgc.annotate([QUERY], pfam), pd.read_parquet(QUERY_TABLE))


def test_bgc_search_returns_linked_spectra(pfam, tmp_path):
    from mibigxfinder.config import BgcSearchConfig
    from mibigxfinder.services import search_bgcs

    cfg = BgcSearchConfig(
        query=[QUERY], reference=REFERENCE, bgc_links=FIXTURES / "bgc_links.csv",
        pfam=pfam, top_n=2, out_dir=tmp_path,
    )
    with pytest.warns(UserWarning, match="linked BGCs are not in"):   # BGC0000005, BGC0003163, ...
        df = search_bgcs(cfg)

    assert df["bgc_id"].tolist() == ["BGC0000664", "BGC0000649"]
    assert df["spectrum_idx"].tolist() == ["CCMSLIB00000079400", "CCMSLIB00000079401"]
    assert df["mibig_rank"].tolist() == [1, 2]
    assert (tmp_path / "bgc_matches.csv").exists()

    # annotating the query live gives BiG-SCAPE's own distances
    merged = df.merge(pd.read_csv(BIGSCAPE, sep="\t"), on="bgc_id", suffixes=("", "_bigscape"))
    assert len(merged) == len(df)
    assert (merged["distance"] - merged["distance_bigscape"]).abs().max() < 1e-12
