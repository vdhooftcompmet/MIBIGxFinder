# src/mibig_finder/config.py
"""Typed configs for each subcommand. CLI and future FastAPI both build these."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field, field_validator

_SPECTRA_EXTS = {".mgf", ".parquet"}

class SpectraSearchConfig(BaseModel):
    query: Path
    bgc_links: Path  
    database: Path 
    tolerance: float = 0.02                                                      
    similarity: Literal["modified_cosine", "cosine_flash"] = "modified_cosine"
    top_n: int = 2
    sim_threshold: float = 0.7
    clean: bool = False
    out_dir: Path = Path("./out")

    @field_validator("query", "database")
    @classmethod
    def _spectra_file_ok(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"file not found: {v}")
        if v.suffix.lower() not in _SPECTRA_EXTS:
            raise ValueError(f"expected .mgf or .parquet, got {v.suffix}")
        return v

    @field_validator("tolerance")
    @classmethod
    def _tolerance_range(cls, v: float) -> float:
        if not (0 < v <= 1.0):
            raise ValueError("tolerance must be in (0, 1.0]")
        return v

    @field_validator("sim_threshold")
    @classmethod
    def _threshold_range(cls, v: float) -> float:
        if not (0 <= v <= 1):
            raise ValueError("sim_threshold must be in [0, 1]")
        return v

    @field_validator("bgc_links")
    @classmethod
    def _file_ok(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"file not found: {v}")
        if v.suffix.lower() not in {".csv", ".tsv"}:
            raise ValueError(f"expected .csv or .tsv, got {v.suffix}")
        return v


class BgcSearchConfig(BaseModel):
    query: list[Path]                                # region .gbk files and/or folders of them
    reference: Path = Path("data/bgc_reference/mibig_5.0.parquet")
    bgc_links: Path = Path("data/library/bgc_links.csv")
    pfam: Path = Path("data/pfam/Pfam-A.hmm")
    top_n: int = 5                                   # closest MIBiG BGCs with spectra, per query
    max_distance: float = 1.0                        # 1.0 = no shared domains
    cores: int = 1
    out_dir: Path = Path("./out")

    @field_validator("query")
    @classmethod
    def _gbk_files(cls, v: list[Path]) -> list[Path]:
        files = []
        for p in v:
            if p.is_dir():
                files += sorted(p.glob("*.gbk"))
            elif p.is_file() and p.suffix.lower() == ".gbk":
                files.append(p)
            else:
                raise ValueError(f"not a .gbk file or folder: {p}")
        if not files:
            raise ValueError("no .gbk files found")
        stems = [f.stem for f in files]
        dupes = sorted({s for s in stems if stems.count(s) > 1})
        if dupes:
            raise ValueError(f"query file names must be unique: {dupes}")
        return files

    @field_validator("reference", "bgc_links", "pfam")
    @classmethod
    def _exists(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"file not found: {v}")
        return v

    @field_validator("max_distance")
    @classmethod
    def _distance_range(cls, v: float) -> float:
        if not (0 < v <= 1):
            raise ValueError("max_distance must be in (0, 1]")
        return v


class StructureSearchConfig(BaseModel):
    query: Path
    database: Path                                   # local MGF (SMILES pulled from spectrum metadata)
    fps: list[Literal["morgan", "map4", "rdkit"]] = Field(default_factory=lambda: ["morgan"])
    similarity_threshold: float = 0.7
    out_dir: Path = Path("./out")

    @field_validator("query")
    @classmethod
    def _query_ok(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"file not found: {v}")
        if v.suffix.lower() not in {".tsv", ".csv", ".smi", ".txt"}:
            raise ValueError(f"expected .tsv/.csv/.smi/.txt, got {v.suffix}")
        return v

    @field_validator("database")
    @classmethod
    def _db_ok(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"file not found: {v}")
        if v.suffix.lower() != ".mgf":
            raise ValueError(f"expected .mgf, got {v.suffix}")
        return v

    @field_validator("similarity_threshold")
    @classmethod
    def _threshold_range(cls, v: float) -> float:
        if not (0 <= v <= 1):
            raise ValueError("similarity_threshold must be in [0, 1]")
        return v