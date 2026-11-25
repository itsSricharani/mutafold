import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

merged = pd.read_csv("data/processed/P01116_merged.csv")
sub = merged[merged["label"].isin(["pathogenic", "benign"])]
p = sub[sub.label == "pathogenic"]
b = sub[sub.label == "benign"]

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for a, col in zip(ax, ["rsa", "plddt"]):
    a.boxplot([b[col], p[col]], tick_labels=["benign", "pathogenic"])
    a.set_title(col)
plt.tight_layout()
plt.savefig("results/figures/pathogenic_vs_benign.png", dpi=200)
print("saved results/figures/pathogenic_vs_benign.png")