"""Audit saved model probabilities and plot the actual designed sequences."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples/proteinmpnn"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=OUT)
    args = parser.parse_args()
    out = args.run.resolve()
    rows = list(csv.DictReader((out / "candidate_summary.csv").open()))
    saved = np.load(out / "designs/probs/1YCR.npz")
    mask = saved["mask"][0]
    selected = np.take_along_axis(saved["log_probs"], saved["S"][:, :, None], axis=2)[:, :, 0]
    scores = -(selected * mask).sum(axis=1) / mask.sum()
    assert mask.sum() == 10
    assert saved["chain_order"].tolist() == [["B", "A"]]
    for score, row in zip(scores, rows):
        assert abs(float(score) - float(row["score"])) < 0.0001
    checks = {
        "probability_shape": list(saved["probs"].shape),
        "chain_order": saved["chain_order"].tolist(),
        "score_mask_nonzero_zero_based": np.flatnonzero(mask).tolist(),
        "rescored_negative_log_likelihood": scores.tolist(),
        "max_header_rounding_error": float(max(abs(float(s) - float(r["score"])) for s, r in zip(scores, rows))),
        "fixed_residues_and_fixed_chain_checked": True,
    }
    parsed = json.loads((out / "parsed.jsonl").read_text().splitlines()[0])
    alphabet = "ACDEFGHIKLMNPQRSTVWYX"
    for seq in saved["S"]:
        decoded = "".join(alphabet[int(x)] for x in seq)
        assert decoded[13:] == parsed["seq_chain_A"]
        assert all(decoded[i] == parsed["seq_chain_B"][i] for i in [2, 6, 9])
    (out / "probability_audit.json").write_text(json.dumps(checks, indent=2) + "\n")
    sequences = [parsed["seq_chain_B"]] + [r["sequence"] for r in rows]
    labels = ["Native"] + [r["candidate"] for r in rows]
    fig, ax = plt.subplots(figsize=(11, 4.8), constrained_layout=True)
    ax.set(xlim=(0.5, 13.5), ylim=(6.5, -0.5), xlabel="Parsed position in chain B")
    ax.set_xticks(range(1, 14)); ax.set_yticks(range(7), labels)
    for pos in [3, 7, 10]:
        ax.axvspan(pos - 0.5, pos + 0.5, color="#d4e8db", zorder=0)
    for i, seq in enumerate(sequences):
        for j, aa in enumerate(seq):
            color = "#19627e" if aa == sequences[0][j] else "#b35336"
            ax.text(j + 1, i, aa, ha="center", va="center", color=color, fontfamily="monospace", fontsize=15)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title("Actual ProteinMPNN samples: green columns are fixed hotspots")
    fig.savefig(out / "sequence_comparison.png", dpi=180)
    print(json.dumps(checks))


if __name__ == "__main__":
    main()
