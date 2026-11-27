## Install

```bash
git clone https://github.com/rtlortega/MIBIGxFinder.git
cd MIBIGxFinder
conda env create -f environment.yml
conda activate mibigxfinder
```

### Pfam database

BGC search also needs the Pfam database. Download it into `data/pfam`:

```bash
mkdir -p data/pfam
wget -P data/pfam https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip data/pfam/Pfam-A.hmm.gz
```

The first BGC annotation presses the hmm profiles (`Pfam-A.hmm.h3*` files are written), so the folder must be writable.

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

### BGC search

Compare query BGCs against MIBiG with BiG-SCAPE's distance and return the spectra linked to the closest MIBiG BGCs (needs the Pfam database, see Install).

```bash
mibigxfinder bgc \
    tests/fixtures/NC_003888.3.region002.gbk \
    --reference tests/fixtures/bgc_reference_small.parquet \
    --links tests/fixtures/bgc_links.csv \
    --pfam tests/fixtures/pfam_small.hmm \
    --out-dir out/
```

The test reference was built with `tests/fixtures/pfam_small.hmm`, which holds only the 21 Pfam families found in the test region, so the example runs in seconds. For real queries, leave out `--pfam` to use `data/pfam/Pfam-A.hmm`.

Positional arguments:
| arg | description |
|---|---|
| `query` | One or more antiSMASH region `.gbk` files, or folders of them such as an antiSMASH output folder; file names must be unique |

From a folder, only `.gbk` files with `region` or `cluster` in the name (and not `final`) are used so the whole-genome `.gbk` in an antiSMASH output folder is skipped. Files given one by one are always used.

Options:
| flag | default | description |
|---|---|---|
| `--reference` | `data/bgc_reference/mibig_5.0.parquet` | BGC reference table (see `build-reference`) |
| `--links` | `data/library/bgc_links.csv` | BGC-links CSV/TSV (see below) |
| `--pfam` | `data/pfam/Pfam-A.hmm` | Pfam file the reference was built with; checked against the reference's `.json` manifest |
| `--top-n` | `5` | Closest MIBiG BGCs with spectra kept per query |
| `--max-distance` | `1.0` | Keep BGCs with a BiG-SCAPE distance below this; `1.0` only drops BGCs without shared domains (BiG-SCAPE groups BGCs into families at `0.3`) |
| `--cores` | `1` | CPU cores for the domain search |
| `--out-dir` | `./out` | Where to write `bgc_matches` and `bgc_summary` (`.csv` / `.parquet`) |

Output columns: `query_id, bgc_id, distance, jaccard, adjacency, dss, mibig_rank`, plus the columns of the links file, one row per query, MIBiG BGC and linked spectrum. `mibig_rank` is the rank among all MIBiG BGCs, so rank 3 means two closer MIBiG BGCs have no spectra.

The terminal only shows a summary, which is also saved as `bgc_summary.csv` / `.parquet`: one line per query with its closest MIBiG BGC with spectra (`closest_bgc`, `distance`, `mibig_rank`, and `compound` if the links file has `bgc_compound`), and how many MIBiG BGCs (`n_bgcs`) and spectra (`n_spectra`) it matched.

### Build the BGC reference

Annotate all MIBiG BGCs with Pfam once; BGC search compares query BGCs against this table (needs the Pfam database, see Install).

```bash
mibigxfinder build-reference --mibig-version 5.0 --cores 4
```

Without `--gbk-dir`, the antiSMASH-processed MIBiG release is downloaded from https://dl.secondarymetabolites.org/mibig/ into `data/bgc_reference/MIBiG/` (MIBiG 5.0 is not available there yet; until then use `--mibig-version 4.0` or pass your own files with `--gbk-dir`). Writes `data/bgc_reference/mibig_5.0.parquet` (one row per gene and Pfam domain) and `mibig_5.0.json` (MIBiG version, Pfam checksum, BiG-SCAPE version and settings). Expect roughly an hour on one core for MIBiG 4.0's ~2,600 BGCs.

### BGC-links file format

Plain CSV, at least these two columns (other columns are ignored):
```csv
spectrum_idx,bgc_id,tanimoto
7,BGC0001082,1.0
7,BGC0003163,0.5
11,BGC0002033,0.7
```
`spectrum_idx` must match the library spectra's `SPECTRUM_IDX` (or `SPECTRUM_ID`) metadata.