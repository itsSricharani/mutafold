"""py3Dmol helpers: color by pLDDT and mark variants."""
import py3Dmol
from Bio.PDB import PDBParser
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
from .fetch import RAW

# Standard AlphaFold confidence colors
BANDS = [(90, 101, "#0053D6"), (70, 90, "#65CBF3"),
         (50, 70, "#FFDB13"), (0, 50, "#FF7D45")]
LABEL_COLORS = {"pathogenic": "red", "benign": "limegreen",
                "uncertain": "orange", "conflicting": "purple"}


def load_structure(acc: str):
    return PDBParser(QUIET=True).get_structure(acc, RAW / f"{acc}_alphafold.pdb")


def residue_info(acc: str) -> dict:
    """{residue number: (one-letter code, pLDDT)} from the AlphaFold file."""
    info = {}
    for res in load_structure(acc).get_residues():
        aa = THREE_TO_ONE.get(res.get_resname().upper(), "X")
        info[res.id[1]] = (aa, res["CA"].get_bfactor())  # B-factor = pLDDT
    return info


def show(acc: str, variants=None, labels=("pathogenic",), color_by="plddt",
         width=800, height=600):
    pdb = (RAW / f"{acc}_alphafold.pdb").read_text()
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb, "pdb")

    if color_by == "plddt":
        info = residue_info(acc)
        for lo, hi, color in BANDS:
            resi = [r for r, (_, p) in info.items() if lo <= p < hi]
            if resi:
                view.setStyle({"resi": resi}, {"cartoon": {"color": color}})
    else:
        view.setStyle({"cartoon": {"color": "lightgrey"}})

    if variants is not None:
        for label in labels:
            sub = variants[variants["label"] == label]
            for pos in sub["position"].unique():
                view.addStyle({"resi": int(pos), "atom": "CA"},
                              {"sphere": {"radius": 1.0,
                                          "color": LABEL_COLORS[label]}})
    view.zoomTo()
    return view


CLUSTER_COLORS = {-1: "#999999", 0: "red", 1: "blue", 2: "limegreen", 3: "orange"}


def show_clusters(acc: str, clusters, width=800, height=600):
    """Color the pathogenic positions by cluster. clusters = DataFrame(position, cluster)."""
    pdb = (RAW / f"{acc}_alphafold.pdb").read_text()
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb, "pdb")
    view.setStyle({"cartoon": {"color": "lightgrey", "opacity": 0.6}})
    for _, row in clusters.iterrows():
        color = CLUSTER_COLORS.get(int(row["cluster"]), "black")
        view.addStyle({"resi": int(row["position"]), "atom": "CA"},
                      {"sphere": {"radius": 1.3, "color": color}})
        view.addLabel(str(int(row["position"])),
                      {"fontSize": 11, "fontColor": "black", "backgroundOpacity": 0.4},
                      {"resi": int(row["position"]), "atom": "CA"})
    view.zoomTo()
    return view