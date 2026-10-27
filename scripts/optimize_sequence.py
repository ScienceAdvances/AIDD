"""Audit a constrained p53 mutation library with official ProteinMPNN weights.

This is fixed-backbone model scoring, not an affinity or expression assay.
"""
from pathlib import Path
import argparse
import copy
import csv
import hashlib
import json
import platform
import random
import sys
import time

import numpy as np
import torch
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT / "data/aidd/ProteinMPNN"
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
CANONICAL = ALPHABET[:-1]
NATIVE = "ETFSDLWKLLPEN"
HOTSPOTS = (3, 7, 10)
SEED = 20261016
COMMIT = "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57"
DEVICE = "cpu"
DESIGN_RUN = ROOT / "examples/proteinmpnn"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mutations(sequence, parent=NATIVE):
    if len(sequence) != len(parent) or any(aa not in CANONICAL for aa in sequence):
        raise ValueError("Expected a canonical 13-residue sequence")
    return [f"{old}{i + 1}{new}" for i, (old, new) in enumerate(zip(parent, sequence)) if old != new]


def properties(sequence):
    changes = mutations(sequence)
    reasons = []
    if any(sequence[i - 1] != NATIVE[i - 1] for i in HOTSPOTS):
        reasons.append("hotspot_changed")
    if len(changes) > 2:
        reasons.append("mutation_budget_exceeded")
    if "C" in sequence:
        reasons.append("new_cysteine")
    if any(sequence[i - 1] == "P" and NATIVE[i - 1] != "P" for i in range(4, 10)):
        reasons.append("new_proline_in_review_region")
    analysis = ProteinAnalysis(sequence)
    return {"mutation_count": len(changes), "mutations_parsed": ";".join(changes),
            "mutations_pdb": ";".join(f"{NATIVE[i]}{i + 17}{aa}" for i, aa in enumerate(sequence) if aa != NATIVE[i]),
            "gravy": analysis.gravy(), "charge_pH7_free_termini": analysis.charge_at_pH(7.0),
            "molecular_weight_free_termini": analysis.molecular_weight(),
            "new_cysteines": sequence.count("C"),
            "new_proline_review_positions": ";".join(str(i) for i in range(4, 10)
                if sequence[i - 1] == "P" and NATIVE[i - 1] != "P"),
            "hard_pass": not reasons, "exclusion_reasons": ";".join(reasons)}


def library():
    rows = [{"candidate_id": "WT", "source": "native", "sequence": NATIVE}]
    for pos in range(1, 14):
        if pos in HOTSPOTS:
            continue
        for aa in CANONICAL:
            if aa == NATIVE[pos - 1]:
                continue
            sequence = NATIVE[:pos - 1] + aa + NATIVE[pos:]
            rows.append({"candidate_id": f"single_{NATIVE[pos - 1]}{pos}{aa}",
                         "source": "single_mutant", "sequence": sequence})
    previous = DESIGN_RUN / "candidate_summary.csv"
    with previous.open(newline="") as stream:
        for row in csv.DictReader(stream):
            rows.append({"candidate_id": "mpnn_" + row["candidate"],
                         "source": "debian_mpnn_sample", "sequence": row["sequence"]})
    if len({r['sequence'] for r in rows}) != len(rows):
        raise ValueError("Duplicate library sequence")
    for row in rows:
        row.update(properties(row["sequence"]))
    return rows


def pareto_flags(rows):
    """Minimize NLL, GRAVY and mutation burden among hard-pass candidates."""
    fields = ["mpnn_nll_mean", "gravy", "mutation_count"]
    for row in rows:
        row["pareto"] = bool(row["hard_pass"]) and not any(
            other["hard_pass"] and
            all(other[key] <= row[key] for key in fields) and
            any(other[key] < row[key] for key in fields)
            for other in rows)


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def features():
    sys.path.insert(0, str(REPO))
    from protein_mpnn_utils import tied_featurize
    parsed = json.loads((DESIGN_RUN / "parsed.jsonl").read_text().splitlines()[0])
    if parsed["seq_chain_B"] != NATIVE or len(parsed["seq_chain_A"]) != 85:
        raise ValueError("Unexpected structure sequence")
    chain_dict = {parsed["name"]: [["B"], ["A"]]}
    fixed = {parsed["name"]: {"B": list(HOTSPOTS), "A": []}}
    tensors = tied_featurize([copy.deepcopy(parsed)], torch.device("cpu"), chain_dict, fixed)
    X, S, mask = tensors[:3]
    chain_M, chain_encoding = tensors[4:6]
    chain_order, chain_M_pos, residue_idx = tensors[6], tensors[10], tensors[12]
    design_mask = mask * chain_M * chain_M_pos
    expected = [i for i in range(13) if i + 1 not in HOTSPOTS]
    if chain_order != [["B", "A"]] or torch.where(design_mask[0] > 0)[0].tolist() != expected:
        raise ValueError("Chain order or design mask changed")
    decoded = "".join(ALPHABET[i] for i in S[0].tolist())
    if decoded != NATIVE + parsed["seq_chain_A"] or int(mask.sum()) != 98:
        raise ValueError("Sequence mapping or coordinate validity failed")
    return X, S, mask, design_mask, residue_idx, chain_encoding


def score_library(rows, repeats, batch_size):
    from protein_mpnn_utils import ProteinMPNN
    checkpoint_path = REPO / "vanilla_model_weights/v_48_020.pt"
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    model = ProteinMPNN(ca_only=False, num_letters=21, node_features=128,
        edge_features=128, hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
        augment_eps=0.0, k_neighbors=checkpoint["num_edges"])
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(DEVICE)
    model.eval()
    X, native_S, mask, design_mask, residue_idx, chain_encoding = features()
    X, native_S, mask, design_mask, residue_idx, chain_encoding = [
        tensor.to(DEVICE) for tensor in (X, native_S, mask, design_mask, residue_idx, chain_encoding)]
    sequences = native_S.repeat(len(rows), 1)
    for index, row in enumerate(rows):
        sequences[index, :13] = torch.tensor([ALPHABET.index(aa) for aa in row["sequence"]], device=DEVICE)
    if not torch.all(sequences[:, 13:] == native_S[:, 13:]):
        raise ValueError("Receptor chain changed")
    for pos in HOTSPOTS:
        if not torch.all(sequences[:, pos - 1] == native_S[:, pos - 1]):
            raise ValueError("Hotspot changed")
    # Identical random draws give the same decoding orders for every candidate.
    generator = torch.Generator().manual_seed(SEED)
    random_orders = torch.randn((repeats, 98), generator=generator)
    nll = np.empty((len(rows), repeats, 98), dtype=np.float32)
    with torch.inference_mode():
        for start in range(0, len(rows), batch_size):
            stop = min(start + batch_size, len(rows))
            count = stop - start
            expand = lambda value: value.repeat(count, *([1] * (value.ndim - 1)))
            for repeat in range(repeats):
                log_probs = model(expand(X), sequences[start:stop], expand(mask),
                    expand(design_mask), expand(residue_idx), expand(chain_encoding),
                    random_orders[repeat:repeat + 1].repeat(count, 1).to(DEVICE))
                selected = torch.gather(log_probs, 2, sequences[start:stop, :, None]).squeeze(-1)
                nll[start:stop, repeat] = -selected.cpu().numpy()
            print(f"Scored candidates {start + 1}-{stop} of {len(rows)}", flush=True)
    score_mask = design_mask[0].cpu().numpy()
    scores = (nll * score_mask).sum(axis=2) / score_mask.sum()
    if not np.isfinite(scores).all():
        raise ValueError("Non-finite model scores")
    return scores, nll, sequences.cpu().numpy(), score_mask, random_orders.numpy()


def plot(rows, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    groups = [("Excluded", [r for r in rows if not r["hard_pass"]], "#b9bdc5"),
              ("Eligible", [r for r in rows if r["hard_pass"] and not r["pareto"]], "#2d7191"),
              ("Pareto set", [r for r in rows if r["pareto"]], "#c06837")]
    for label, group, color in groups:
        ax.scatter([r["gravy"] for r in group], [r["mpnn_nll_mean"] for r in group],
                   s=32, alpha=0.8, c=color, label=label)
    wt = rows[0]
    ax.scatter([wt["gravy"]], [wt["mpnn_nll_mean"]], marker="*", s=170, c="#172b44", label="WT")
    for row in rows:
        if row.get("shortlisted"):
            ax.annotate(row["candidate_id"], (row["gravy"], row["mpnn_nll_mean"]),
                        xytext=(4, 5), textcoords="offset points", fontsize=8)
    ax.set(xlabel="GRAVY (sequence hydropathy proxy)", ylabel="Mean ProteinMPNN NLL (10 positions)",
           title="Actual fixed-backbone scores: candidate selection is a tradeoff")
    ax.legend(fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "tradeoffs.png", dpi=180)
    plt.close(fig)


def main():
    global DEVICE, DESIGN_RUN
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-run", type=Path, default=DESIGN_RUN)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    args = parser.parse_args()
    DEVICE = args.device
    DESIGN_RUN = args.design_run.resolve()
    if DEVICE == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    if args.repeats < 2 or args.batch_size < 1:
        parser.error("Use at least two repeats and a positive batch size")
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    output.mkdir(parents=True)
    started = time.monotonic()
    torch.set_num_threads(4)
    torch.manual_seed(SEED)
    random.seed(SEED)
    np.random.seed(SEED)
    # Features also add the frozen official source directory to the import path.
    features()
    rows = library()
    with (output / "library.fasta").open("w") as stream:
        for row in rows:
            stream.write(f">{row['candidate_id']}\n{row['sequence']}\n")
    scores, nll, sequences, score_mask, random_orders = score_library(rows, args.repeats, args.batch_size)
    for index, row in enumerate(rows):
        paired = scores[index] - scores[0]
        row.update(mpnn_nll_mean=float(scores[index].mean()),
                   decoding_order_sd=float(scores[index].std(ddof=1)),
                   delta_nll_vs_wt=float(paired.mean()),
                   delta_nll_order_sd=float(paired.std(ddof=1)),
                   fraction_orders_better_than_wt=float(np.mean(paired < 0)))
    pareto_flags(rows)
    # Keep WT and three eligible Pareto candidates separated by >=2 substitutions.
    shortlist = [rows[0]]
    for row in sorted((r for r in rows if r["pareto"] and r["candidate_id"] != "WT"),
                      key=lambda r: (r["mpnn_nll_mean"], r["candidate_id"])):
        if all(sum(a != b for a, b in zip(row["sequence"], other["sequence"])) >= 2
               for other in shortlist[1:]):
            shortlist.append(row)
        if len(shortlist) == 4:
            break
    for row in rows:
        row["shortlisted"] = row in shortlist
    write_csv(output / "candidates.csv", rows)
    write_csv(output / "shortlist.csv", shortlist)
    np.savez_compressed(output / "score_audit.npz", scores=scores, per_position_nll=nll,
                        S=sequences, score_mask=score_mask, random_orders=random_orders,
                        candidate_ids=np.array([r["candidate_id"] for r in rows]))
    imports = {name: __import__(name).__version__ for name in ["numpy", "torch", "Bio", "matplotlib"]}
    paths = [ROOT / "data/peptide/1YCR_original.pdb",
             DESIGN_RUN / "parsed.jsonl",
             DESIGN_RUN / "candidate_summary.csv",
             REPO / "protein_mpnn_utils.py", REPO / "vanilla_model_weights/v_48_020.pt",
             Path(__file__)]
    summary = {"completed": True, "design_run": str(DESIGN_RUN), "date": __import__("datetime").date.today().isoformat(), "device": DEVICE,
        "platform": platform.platform(), "python": platform.python_version(), "packages": imports,
        "official_commit": COMMIT, "model": "v_48_020", "seed": SEED,
        "decoding_orders": args.repeats, "common_orders_for_all_candidates": True,
        "chain_order": ["B", "A"], "score_mask_positions_B": [int(i + 1) for i in np.flatnonzero(score_mask)],
        "fixed_hotspots_B": list(HOTSPOTS), "receptor_and_hotspots_checked": True,
        "score_units": "dimensionless average negative natural log probability",
        "scoring_temperature": "not applied; raw model log probabilities",
        "candidate_count": len(rows), "single_mutants": sum(r["source"] == "single_mutant" for r in rows),
        "eligible_count": sum(r["hard_pass"] for r in rows), "pareto_count": sum(r["pareto"] for r in rows),
        "shortlist": [r["candidate_id"] for r in shortlist],
        "hard_rules": {"max_mutations": 2, "no_new_cysteine": True, "no_new_proline_positions": list(range(4, 10))},
        "pareto_objectives_minimize": ["mpnn_nll_mean", "gravy", "mutation_count"],
        "input_and_script_sha256": {str(p.relative_to(ROOT)): digest(p) for p in paths},
        "source_manifest_sha256": digest(ROOT / "source_manifest.json"),
        "elapsed_seconds": time.monotonic() - started,
        "structural_reprediction_completed": False, "experimental_validation": False}
    plot(rows, output)
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
