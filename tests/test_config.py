"""Which query files BGC search takes from its arguments (no BiG-SCAPE needed)."""

import pytest
from pydantic import ValidationError

from mibigxfinder.config import BgcSearchConfig


def _config(tmp_path, query):
    for name in ("reference.parquet", "links.tsv", "Pfam-A.hmm"):   # only checked to exist
        (tmp_path / name).touch()
    return BgcSearchConfig(
        query=query, reference=tmp_path / "reference.parquet",
        bgc_links=tmp_path / "links.tsv", pfam=tmp_path / "Pfam-A.hmm",
    )


def test_antismash_folder_gives_only_region_files(tmp_path):
    run = tmp_path / "antismash_run"
    run.mkdir()
    for name in [
        "NC_1.region001.gbk", "NC_1.region002.gbk",   # antiSMASH 5+ regions
        "NC_1.cluster001.gbk",                        # antiSMASH 4 cluster
        "NC_1.gbk", "NC_1.final.gbk",                 # whole genome (antiSMASH 5+ / 4)
        "NC_1.region003.json",
    ]:
        (run / name).touch()

    cfg = _config(tmp_path, [run])
    assert [f.name for f in cfg.query] == ["NC_1.cluster001.gbk", "NC_1.region001.gbk", "NC_1.region002.gbk"]


def test_files_given_one_by_one_are_always_used(tmp_path):
    gbk = tmp_path / "my_bgc.gbk"
    gbk.touch()
    assert _config(tmp_path, [gbk]).query == [gbk]


def test_folder_without_region_files_is_an_error(tmp_path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "NC_1.gbk").touch()
    with pytest.raises(ValidationError, match="no region .gbk files"):
        _config(tmp_path, [genomes])
