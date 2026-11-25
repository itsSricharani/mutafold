"""Per-residue structural features: secondary structure, burial, pLDDT."""
import pandas as pd
import pydssp
from Bio.PDB.SASA import ShrakeRupley
from .fetch import RAW
from .viz import load_structure, THREE_TO_ONE

# Theoretical maximum solvent accessibility per residue type (Tien et al. 2013)
MAX_ASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223,
           "Q": 225, "G": 104, "H": 224, "I": 197, "L": 201, "K": 236,
           "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
           "Y": 263, "V": 174}


def per_residue_table(acc: str) -> pd.DataFrame:
    structure = load_structure(acc)
    residues = list(structure.get_residues())

    # Solvent accessibility (Shrake-Rupley), per residue
    ShrakeRupley().compute(structure, level="R")

    # Secondary structure: H = helix, E = strand, - = loop
    coord = pydssp.read_pdbtext((RAW / f"{acc}_alphafold.pdb").read_text())
    ss = pydssp.assign(coord, out_type="c3")

    rows = []
    for res, s in zip(residues, ss):
        aa = THREE_TO_ONE.get(res.get_resname().upper(), "X")
        rows.append({
            "position": res.id[1],
            "aa": aa,
            "plddt": res["CA"].get_bfactor(),
            "ss": s,
            "rsa": min(res.sasa / MAX_ASA[aa], 1.0) if aa in MAX_ASA else None,
        })
    return pd.DataFrame(rows)