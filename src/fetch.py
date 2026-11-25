"""Fetch structure, variants, and annotations for a UniProt protein."""
from pathlib import Path
import json
import requests
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
VARIATION_API = "https://www.ebi.ac.uk/proteins/api/variation/{acc}"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{acc}.json"


def _get(url, **kwargs):
    r = requests.get(url, timeout=60, **kwargs)
    r.raise_for_status()
    return r


def fetch_alphafold(acc: str) -> Path:
    """Download the AlphaFold PDB for a UniProt accession."""
    RAW.mkdir(parents=True, exist_ok=True)
    meta = _get(ALPHAFOLD_API.format(acc=acc)).json()[0]
    (RAW / f"{acc}_alphafold_meta.json").write_text(json.dumps(meta, indent=2))
    pdb_path = RAW / f"{acc}_alphafold.pdb"
    pdb_path.write_bytes(_get(meta["pdbUrl"]).content)
    print(f"AlphaFold: saved {pdb_path.name} (entry {meta.get('entryId')})")
    return pdb_path


def fetch_variants(acc: str) -> Path:
    """Download all variants (with ClinVar-derived clinical significance)."""
    RAW.mkdir(parents=True, exist_ok=True)
    data = _get(VARIATION_API.format(acc=acc),
                headers={"Accept": "application/json"}).json()
    path = RAW / f"{acc}_variants.json"
    path.write_text(json.dumps(data))
    print(f"Variants: saved {path.name} ({len(data.get('features', []))} features)")
    return path


def fetch_uniprot(acc: str) -> Path:
    """Download the UniProt entry (domains, binding sites, etc.)."""
    RAW.mkdir(parents=True, exist_ok=True)
    data = _get(UNIPROT_API.format(acc=acc)).json()
    path = RAW / f"{acc}_uniprot.json"
    path.write_text(json.dumps(data))
    print(f"UniProt: saved {path.name}")
    return path


def variants_to_table(acc: str) -> pd.DataFrame:
    """Flatten the variant JSON into a tidy missense table."""
    data = json.loads((RAW / f"{acc}_variants.json").read_text())
    rows = []
    for f in data.get("features", []):
        if f.get("type") != "VARIANT":
            continue
        sig = f.get("clinicalSignificances") or []
        sig_types = sorted({
            (s.get("type") if isinstance(s, dict) else str(s)).lower()
            for s in sig if s
        })
        rows.append({
            "position": int(f["begin"]) if f.get("begin") else None,
            "wild_type": f.get("wildType"),
            "mutant": f.get("alternativeSequence"),
            "consequence": f.get("consequenceType"),
            "clinical_significance": "; ".join(sig_types),
            "n_sources": len(f.get("xrefs") or []),
        })
    df = pd.DataFrame(rows).dropna(subset=["position"])
    df["position"] = df["position"].astype(int)
    return df.sort_values("position").reset_index(drop=True)


def fetch_all(acc: str):
    fetch_alphafold(acc)
    fetch_variants(acc)
    fetch_uniprot(acc)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    df = variants_to_table(acc)
    out = PROCESSED / f"{acc}_variants.csv"
    df.to_csv(out, index=False)
    print(f"Table: saved {out.name} ({len(df)} variants)")
    return df