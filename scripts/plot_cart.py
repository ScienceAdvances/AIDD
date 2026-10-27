"""Plot original observations and descriptive mean/SEM; do not invent donor identities."""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples/debian13_20261023/cart"
COLORS = {"CAR[A]": "#248e69", "CAR[B]": "#ba4949", "CAR[C]": "#8e60a5", "CAR[D]": "#be7b36", "CAR[E]": "#2868b1", "Mock": "#6f7580"}


def panel(ax, data, name, groups, title, ylabel):
    sub = data[data.panel == name]
    for i, g in enumerate(groups):
        v = sub.loc[sub.group == g, "value"].to_numpy()
        if not len(v):
            ax.text(i, 0.92, "n.a.", transform=ax.get_xaxis_transform(), ha="center", color="#6f7580")
            continue
        c = COLORS.get(g, "#a85d54")
        offsets = np.linspace(-0.15, 0.15, len(v))
        ax.scatter(i + offsets, v, color=c, s=24, alpha=0.85, zorder=3)
        sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0
        ax.errorbar(i, v.mean(), yerr=sem, color="black", capsize=4, fmt="_", ms=18, zorder=4)
    ax.set_xticks(range(len(groups)), groups, rotation=25, ha="right")
    ax.set(title=title, ylabel=ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axhline(0, color="#c0c6ca", lw=0.7)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=OUT / "source_observations.csv")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(args.input)
    groups = ["CAR[A]", "CAR[B]", "CAR[C]", "CAR[D]", "CAR[E]", "Mock"]
    fig, a = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    panel(a[0], d, "1E", groups, "Fig. 1E: U-87MG", "Relative target viability (%)")
    panel(a[1], d, "1F", groups, "Fig. 1F: U-343MG", "Relative target viability (%)")
    fig.suptitle("Same CAR set, different target-cell backgrounds\nSource observations with mean +/- SEM", fontsize=12)
    fig.savefig(out / "target_background_comparison.png", dpi=180); plt.close(fig)
    fig, a = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    panel(a[0], d, "2B", ["CAR[A]", "CAR[B]", "CAR[E]"], "Fig. 2B: surface expression", "CAR expression (MFI)")
    panel(a[1], d, "2I", ["CAR[A]", "CAR[B]", "CAR[E]", "Mock"], "Fig. 2I: unstimulated IFN-gamma", "Source workbook values (unit not stated)")
    fig.savefig(out / "expression_and_background.png", dpi=180); plt.close(fig)
    groups = ["CAR[B]", "[B]-H2", "[B]-L3", "[B]-H2L3wt", "CAR[E]", "Mock"]
    fig, a = plt.subplots(1, 3, figsize=(14, 4.7), constrained_layout=True)
    panel(a[0], d, "4H", groups, "Fig. 4H: background signal", "Relative IFN-gamma level")
    panel(a[1], d, "S7B", groups, "Fig. S7B: antigen binding", "Antigen binding (MFI)")
    panel(a[2], d, "S7I_5:1", groups, "Fig. S7I: E:T = 5:1", "Relative target viability (%)")
    fig.suptitle("Reducing background signal must preserve binding and function", fontsize=13)
    fig.savefig(out / "mutation_tradeoffs.png", dpi=180); plt.close(fig)
    print("Saved 3 figures from original source observations")


if __name__ == "__main__":
    main()
