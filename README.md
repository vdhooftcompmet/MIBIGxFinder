## Install

```bash
conda create -n mibigxfinder python=3.12 -y
conda activate mibigxfinder
git clone https://github.com/rtlortega/MIBIGxFinder.git
cd mibigxfinder
pip install -e .
```

### BGC search (optional)

BGC search uses the distance metric of [BiG-SCAPE 2](https://github.com/medema-group/BiG-SCAPE) and needs BiG-SCAPE (2.0.3) and the Pfam database. Install BiG-SCAPE from bioconda (Linux or macOS) into the same environment:

```bash
conda install -n mibigxfinder -c conda-forge -c bioconda bigscape=2.0.3 -y
```

Download Pfam into `data/pfam` (~2 GB):

```bash
mkdir -p data/pfam
wget -P data/pfam https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip data/pfam/Pfam-A.hmm.gz
```

The first BGC annotation presses the profiles (`Pfam-A.hmm.h3*` files are written), so the folder must be writable.

## Usage

### Spectra search

Match query MS/MS spectra against the reference library and annotate hits with BGC IDs.

```bash
mibigxfinder spectra \
    tests/fixtures/query_fake.mgf \
    tests/fixtures/library_fake.mgf \
    tests/fixtures/bgc_links.csv \
    --out-dir out/
```

Positional arguments:
| arg | description |
|---|---|
| `query` | Query spectra — `.mgf` or `.parquet` |
| `database` | Reference library spectra — `.mgf` or `.parquet` |
| `bgc_links` | CSV mapping `spectrum_idx,bgc_id` (one row per library-spectrum / BGC pair; same spectrum may appear multiple times) |

Options:
| flag | default | description |
|---|---|---|
| `--similarity` | `modified_cosine` | `modified_cosine` or `cosine` |
| `--tolerance` | `0.02` | Fragment peak tolerance (Da) |
| `--top-n` | `2` | Top-N matches kept per query, before expansion |
| `--threshold` | `0.7` | Minimum similarity to keep |
| `--out-dir` | `./out` | Where to write `spectra_matches.csv` / `.parquet` |

Output columns: `feature_id, library_id, bgc_id, inchikey, similarity, matched_peaks, rank`. Rank is re-numbered after the BGC join, so a library hit with N BGCs produces N rows ranked 1..N for that query.

### Structure search

```bash
mibigxfinder structure query_smiles.csv tests/fixtures/library_fake.mgf --out-dir out/
```

### BGC-links file format

Plain CSV, at least these two columns (other columns are ignored):
```csv
spectrum_idx,bgc_id,tanimoto
7,BGC0001082,1.0
7,BGC0003163,0.5
11,BGC0002033,0.7
```
`spectrum_idx` must match the library spectra's `SPECTRUM_IDX` (or `SPECTRUM_ID`) metadata.