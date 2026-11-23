"""Runs MIBiG-Finder on MIBiG-Finder/tests/fixtures, to see how the Finder works.

Usage (from anywhere):
    python run_fixtures.py [out_dir]     # default out_dir: ./finder_out

Only shim: maps Scores.to_array('score'|'matches') to the '<Similarity>_score|_matches'
fields of matchms 0.27.0. Without it the Finder crashes (see PROJECT_LOG.md).
"""
import sys, warnings, logging
warnings.filterwarnings("ignore")
logging.getLogger("matchms").setLevel(logging.ERROR)
from pathlib import Path
import pandas as pd
from matchms.Scores import Scores

# shim, a small patch that makes the Finder run on matchms 0.27.0. It only matters for the spectra search. 
_orig = Scores.to_array # keep matchms' original method
def to_array(self, name=None): # define a replacement
    if name in ("score", "matches"):
        name = [n for n in self.score_names if n.endswith("_" + name)][0]
    return _orig(self, name) # then call the original
Scores.to_array = to_array  # swap the replacement in

# purpuse of above lines: If the Finder asks for "score", the replacement looks through the real field names and finds the one ending in _score, which is ModifiedCosine_score.


from mibig_finder.config import SpectraSearchConfig, StructureSearchConfig
from mibig_finder.services import search_spectra, search_structures

HERE = Path(__file__).resolve().parent
fx = HERE #/ "MIBiG-Finder" / "tests" / "fixtures"
out = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "finder_out"
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_colwidth", 40)

print("=== SPECTRA SEARCH: query_fake.mgf vs library_fake.mgf ===")
df = search_spectra(SpectraSearchConfig(query=fx/"query_fake.mgf", database=fx/"library_fake.mgf", out_dir=out/"spectra"))
print(df.to_string(index=False))

print("\n=== STRUCTURE SEARCH: structures_fake.tsv vs library_fake.mgf (morgan, threshold 0.7) ===")
df2 = search_structures(StructureSearchConfig(query=fx/"structures_fake.tsv", database=fx/"library_fake.mgf", out_dir=out/"structure"))
print(df2.to_string(index=False))

print("\n=== STRUCTURE SEARCH: smiles.csv ===")
df3 = search_structures(StructureSearchConfig(query=fx/"smiles.csv", database=fx/"library_fake.mgf", out_dir=out/"structure_csv"))
print(df3.to_string(index=False) if len(df3) else "(no rows)")

print(f"\nOutputs written to {out}")
