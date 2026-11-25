import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from src.variants import load_missense, validate_against_structure
from src.viz import show, show_clusters

ACC = "P01116"
st.set_page_config(page_title="MutaFold", layout="wide")



@st.cache_data
def load_data():
    df = validate_against_structure(load_missense(ACC), ACC)
    clusters = pd.read_csv(ROOT / "data" / "processed" / "P01116_clusters_eps5.5.csv")
    return df, clusters


df, clusters = load_data()

st.title("MutaFold: KRAS variants on the AlphaFold structure")
st.caption("Educational research tool, not for diagnosis. Structure: AlphaFold DB. "
           "Variants: EBI Proteins API (ClinVar-derived).")

# ---- Sidebar ----
st.sidebar.header("Display")
mode = st.sidebar.radio("View", ["Variant classes", "3D clusters (pathogenic)"])

if mode == "Variant classes":
    chosen = st.sidebar.multiselect(
        "Show variant classes",
        ["pathogenic", "benign", "uncertain", "conflicting"],
        default=["pathogenic", "benign"],
    )
    color_by = st.sidebar.radio("Backbone color", ["plddt", "none"],
                                format_func=lambda x: "AlphaFold confidence (pLDDT)"
                                if x == "plddt" else "Plain grey")
    view = show(ACC, df, labels=tuple(chosen), color_by=color_by,
                width=600, height=600)
else:
    view = show_clusters(ACC, clusters, width=600, height=600)


SCROLL_TOP = """
<script>
window.addEventListener('load', function () {
  setTimeout(function () {
    try {
      var d = window.parent.document;
      var m = d.querySelector('[data-testid="stMain"]') || d.querySelector('section.main');
      if (m) m.scrollTo(0, 0);
      window.parent.scrollTo(0, 0);
    } catch (e) {}
  }, 300);
});
</script>
"""

# ---- Layout ----
left, right = st.columns([5, 3])

with left:
    st.iframe(view._make_html() + SCROLL_TOP, height=620)
    if mode == "Variant classes":
        st.caption("Red = pathogenic, green = benign, orange = uncertain, "
                   "purple = conflicting. Backbone: blue = very high confidence, "
                   "cyan = confident, yellow = low, orange = very low.")
    else:
        st.caption("Red = core cluster, blue = cluster 1, green = cluster 2, "
                   "grey = not in any cluster (DBSCAN, eps = 5.5 Å).")

with right:
    st.subheader("Variant counts")
    st.dataframe(df["label"].value_counts().rename("count"))

    st.subheader("Pathogenic hotspots")
    hot = (df[df["label"] == "pathogenic"].groupby("position").size()
           .sort_values(ascending=False).head(10))
    st.bar_chart(hot)

    st.subheader("Look up a position")
    pos = st.number_input("Residue number", min_value=1, max_value=189, value=12)
    hits = df[df["position"] == pos][["name", "label", "clinical_significance"]]
    if hits.empty:
        st.write("No variants recorded at this position.")
    else:
        st.dataframe(hits, hide_index=True)
    c = clusters[clusters["position"] == pos]
    if not c.empty:
        cl = int(c["cluster"].iloc[0])
        st.write("Cluster:", "none (outlier)" if cl == -1 else cl)

st.markdown("---")
st.caption("Caveats: AlphaFold predicts the wild-type protein, not mutants. "
           "ClinVar labels reflect which positions researchers studied.")