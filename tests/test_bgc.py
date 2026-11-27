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
        result = search_bgcs(cfg)

    # similarities: every reference BGC sharing a domain, closest first, with BiG-SCAPE's own
    # distances (the query is annotated live here)
    expected = pd.read_csv(BIGSCAPE, sep="\t")
    expected = expected[expected["distance"] < 1].reset_index(drop=True)
    sims = pd.read_csv(tmp_path / "bgc_similarities.csv")
    assert sims.columns.tolist() == ["query_id", "bgc_id", "distance", "jaccard", "adjacency", "dss", "mibig_rank"]
    assert sims["bgc_id"].tolist() == expected["bgc_id"].tolist()
    assert (sims["distance"] - expected["distance"]).abs().max() < 1e-12
    assert sims["mibig_rank"].tolist() == list(range(1, len(sims) + 1))

    # matches: the spectra of the top_n=2 closest BGCs with spectra
    df = result.matches
    assert df.columns.tolist() == ["query_id", "bgc_id", "distance", "mibig_rank", "spectrum_idx", "tanimoto"]
    assert df["bgc_id"].tolist() == ["BGC0000664", "BGC0000649"]
    assert df["spectrum_idx"].tolist() == ["CCMSLIB00000079400", "CCMSLIB00000079401"]
    assert df["mibig_rank"].tolist() == [1, 2]
    assert (tmp_path / "bgc_matches.csv").exists()

    summary = pd.read_csv(tmp_path / "bgc_summary.csv")
    assert summary[["query_id", "closest_bgc", "n_bgcs", "n_spectra"]].values.tolist() == [
        ["NC_003888.3.region002", "BGC0000664", 2, 2]
    ]


def test_summary_has_one_row_per_query():
    from mibigxfinder.services import summarize_bgc_matches

    matches = pd.DataFrame({
        "query_id": ["q1", "q1", "q1", "q2"],
        "bgc_id": ["BGC2", "BGC1", "BGC1", "BGC3"],
        "distance": [0.4, 0.2, 0.2, 0.7],
        "mibig_rank": [2, 1, 1, 3],
        "spectrum_idx": ["s1", "s2", "s3", "s2"],
        "bgc_compound": ["c", "a", "b", "d"],
    })
    s = summarize_bgc_matches(matches, ["q1", "q2", "q3"])   # q3: no matches

    assert s["query_id"].tolist() == ["q1", "q2", "q3"]
    assert s["closest_bgc"].tolist()[:2] == ["BGC1", "BGC3"]
    assert s["compound"].tolist()[:2] == ["a (+1)", "d"]
    assert s["mibig_rank"].tolist()[:2] == [1, 3] and pd.isna(s["mibig_rank"].iloc[2])
    assert s["n_bgcs"].tolist() == [2, 1, 0]
    assert s["n_spectra"].tolist() == [3, 1, 0]
