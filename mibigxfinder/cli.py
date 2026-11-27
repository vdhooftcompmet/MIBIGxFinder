# │   ├── cli.py                 # Typer: spectra | structure | bgc | db
"""mibigxfinder/cli.py"""
import typer
from pathlib import Path
from typing import Optional

app = typer.Typer(help="MIBiG-Finder: match spectra and structures against a MIBiG-linked library.")


@app.command()
def spectra(
    query: Path = typer.Argument(..., help="Query MGF or Parquet file"),
    database: Path = typer.Argument(..., help="Reference library MGF or Parquet file"),
    bgc_links: Path = typer.Argument(..., help="BGC links CSV (spectrum_idx,bgc_id)"),
    similarity: str = typer.Option("modified_cosine", help="Similarity function: modified_cosine or cosine"),
    tolerance: float = typer.Option(0.02, help="Fragment peak tolerance (Da)"),
    top_n: int = typer.Option(2, help="Top-N matches per query"),
    threshold: float = typer.Option(0.7, help="Minimum similarity score"),
    out_dir: Optional[Path] = typer.Option("./out", help="Output directory"),
):
    """Match query spectra against the reference library."""
    from mibigxfinder.config import SpectraSearchConfig
    from mibigxfinder.services import search_spectra

    cfg = SpectraSearchConfig(
        query=query,
        database=database,
        bgc_links=bgc_links,          
        similarity=similarity,
        tolerance=tolerance,
        top_n=top_n,
        sim_threshold=threshold,
        out_dir=out_dir,
    )
    df = search_spectra(cfg)
    typer.echo(df.to_string(index=False))


@app.command()
def structure(
    query: Path = typer.Argument(..., help="SMILES file (CSV with 'smiles' column, or one per line)"),
    database: Path = typer.Argument(..., help="Reference library MGF file"),
    fps: str = typer.Option("morgan", help="Fingerprints, comma-separated: morgan,map4,rdkit"),
    threshold: float = typer.Option(0.7, help="Minimum Tanimoto similarity"),
    out_dir: Optional[Path] = typer.Option("./out", help="Output directory"),
):
    """Match query structures against SMILES in the reference library."""
    from mibigxfinder.config import StructureSearchConfig
    from mibigxfinder.services import search_structures

    cfg = StructureSearchConfig(
        query=query,
        database=database,
        fps=[f.strip() for f in fps.split(",")],
        similarity_threshold=threshold,
        out_dir=out_dir,
    )
    df = search_structures(cfg)
    typer.echo(df.to_string(index=False))


@app.command()
def bgc(
    query: list[Path] = typer.Argument(..., help="Query antiSMASH region .gbk files, or folders of them such as an antiSMASH output folder"),
    reference: Path = typer.Option(Path("data/bgc_reference/mibig_5.0.parquet"), help="BGC reference table (see build-reference)"),
    bgc_links: Path = typer.Option(Path("data/library/bgc_links.csv"), "--links", help="BGC links CSV/TSV (spectrum_idx,bgc_id,...)"),
    pfam: Path = typer.Option(Path("data/pfam/Pfam-A.hmm"), help="Pfam-A.hmm the reference was built with"),
    top_n: int = typer.Option(5, help="Closest MIBiG BGCs with spectra kept per query"),
    max_distance: float = typer.Option(1.0, help="Keep BGCs with a BiG-SCAPE distance below this (1.0 = no shared domains)"),
    cores: int = typer.Option(1, help="CPU cores for the domain search"),
    out_dir: Optional[Path] = typer.Option("./out", help="Output directory"),
):
    """Find MIBiG BGCs similar to query BGCs and return their linked spectra (needs BiG-SCAPE)."""
    from mibigxfinder.config import BgcSearchConfig
    from mibigxfinder.services import search_bgcs

    cfg = BgcSearchConfig(
        query=query,
        reference=reference,
        bgc_links=bgc_links,
        pfam=pfam,
        top_n=top_n,
        max_distance=max_distance,
        cores=cores,
        out_dir=out_dir,
    )
    df = search_bgcs(cfg)
    typer.echo(df.to_string(index=False))


@app.command("build-reference")
def build_reference(
    pfam: Path = typer.Option(Path("data/pfam/Pfam-A.hmm"), help="Pfam-A.hmm file"),
    mibig_version: str = typer.Option("5.0", help="MIBiG version"),
    gbk_dir: Optional[Path] = typer.Option(None, help="Folder of antiSMASH-processed MIBiG .gbk files (default: download)"),
    out_dir: Path = typer.Option(Path("data/bgc_reference"), help="Output directory"),
    cores: int = typer.Option(1, help="CPU cores for the domain search"),
):
    """Build the MIBiG domain table used by BGC search (needs BiG-SCAPE)."""
    from mibigxfinder.bgc import build_reference as build

    out = build(out_dir, pfam, mibig_version=mibig_version, gbk_dir=gbk_dir, cores=cores)
    typer.echo(f"wrote {out} and {out.with_suffix('.json').name}")


if __name__ == "__main__":
    app()