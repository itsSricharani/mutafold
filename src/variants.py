"""Clean and classify the variant table."""
import pandas as pd
from .fetch import PROCESSED
from .viz import residue_info

AA = set("ACDEFGHIKLMNPQRSTVWY")


def classify(sig: str) -> str:
    s = (sig or "").lower()
    has_path = "pathogenic" in s.replace("likely benign", "")
    has_ben = "benign" in s
    if has_path and has_ben:
        return "conflicting"
    if has_path:
        return "pathogenic"
    if has_ben:
        return "benign"
    if "uncertain" in s:
        return "uncertain"
    return "unlabeled"


def load_missense(acc: str) -> pd.DataFrame:
    df = pd.read_csv(PROCESSED / f"{acc}_variants.csv").fillna("")
    df = df[(df["wild_type"].isin(AA)) & (df["mutant"].isin(AA))].copy()
    df["label"] = df["clinical_significance"].apply(classify)
    df["name"] = df["wild_type"] + df["position"].astype(str) + df["mutant"]
    return df.reset_index(drop=True)

def validate_against_structure(df: pd.DataFrame, acc: str) -> pd.DataFrame:
    """Drop variants whose wild-type residue disagrees with the AlphaFold structure."""
    info = residue_info(acc)
    structure_aa = df["position"].map(lambda p: info.get(p, ("?", 0))[0])
    bad = df["wild_type"] != structure_aa
    if bad.any():
        print(f"Dropping {bad.sum()} variants that don't match the structure:")
        print(df.loc[bad, ["name", "label"]].to_string(index=False))
    return df[~bad].reset_index(drop=True)