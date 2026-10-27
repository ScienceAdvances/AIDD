"""Check library rules, score arithmetic, provenance and batch invariance."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import torch
import optimize_sequence as workflow


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    output = args.run.resolve()
    summary = json.loads((output / "run_summary.json").read_text())
    workflow.DEVICE = summary.get("device", "cpu")
    workflow.DESIGN_RUN = Path(summary.get("design_run", workflow.DESIGN_RUN))
    if not summary["completed"]:
        raise ValueError("Incomplete run")
    rows = list(csv.DictReader((output / "candidates.csv").open()))
    saved = np.load(output / "score_audit.npz")
    assert len(rows) == 197
    assert sum(row["source"] == "single_mutant" for row in rows) == 190
    assert len({row["sequence"] for row in rows}) == 197
    assert np.array_equal(saved["candidate_ids"], [row["candidate_id"] for row in rows])
    for row in rows:
        assert all(row["sequence"][i - 1] == workflow.NATIVE[i - 1] for i in workflow.HOTSPOTS)
        if row["source"] == "single_mutant":
            assert len(workflow.mutations(row["sequence"])) == 1
    for relative, expected in summary["input_and_script_sha256"].items():
        assert workflow.digest(workflow.ROOT / relative) == expected
    assert workflow.digest(workflow.ROOT / "source_manifest.json") == summary["source_manifest_sha256"]
    features = workflow.features()
    baseline = features[1][0].numpy()
    assert np.all(saved["S"][:, 13:] == baseline[13:])
    assert saved["score_mask"].sum() == 10
    assert np.array_equal(saved["score_mask"], features[3][0].numpy())
    recalculated = (saved["per_position_nll"] * saved["score_mask"]).sum(axis=2) / 10
    assert np.allclose(recalculated, saved["scores"], atol=1e-6)
    # Check arithmetic against the author's own score function.
    from protein_mpnn_utils import _scores
    selected_log_probs = torch.zeros((len(rows), 98, 21))
    index = torch.from_numpy(saved["S"][:, :, None]).long()
    selected_log_probs.scatter_(2, index, -torch.from_numpy(saved["per_position_nll"][:, 0, :, None]))
    official_scores = _scores(torch.from_numpy(saved["S"]).long(), selected_log_probs,
                             torch.from_numpy(saved["score_mask"]).repeat(len(rows), 1))
    assert np.allclose(official_scores.numpy(), saved["scores"][:, 0], atol=1e-6)
    assert not workflow.properties("ETASDLWKLLPEN")["hard_pass"]
    assert not workflow.properties("CTFSDLWKLLPEN")["hard_pass"]
    assert not workflow.properties("ETFPDLWKLLPEN")["hard_pass"]
    assert workflow.properties(workflow.NATIVE)["hard_pass"]
    fixture = [
        dict(hard_pass=True, mpnn_nll_mean=2.0, gravy=0.0, mutation_count=0),
        dict(hard_pass=True, mpnn_nll_mean=1.0, gravy=-1.0, mutation_count=1),
        dict(hard_pass=True, mpnn_nll_mean=1.2, gravy=0.0, mutation_count=1),
        dict(hard_pass=False, mpnn_nll_mean=0.5, gravy=-3.0, mutation_count=1)]
    workflow.pareto_flags(fixture)
    assert [row["pareto"] for row in fixture] == [True, True, False, False]
    # The original batch contained 16 sequences; rerun two in a smaller batch.
    torch.set_num_threads(4)
    repeated, _, _, _, orders = workflow.score_library(
        rows[:2], repeats=summary["decoding_orders"], batch_size=2)
    assert np.array_equal(orders, saved["random_orders"])
    maximum_error = float(np.max(np.abs(repeated - saved["scores"][:2])))
    assert maximum_error < 2e-5
    verification = {
        "date": __import__("datetime").date.today().isoformat(), "passed": True,
        "candidate_count_and_complete_single_mutation_library": True,
        "input_and_script_hashes_match": True,
        "receptor_and_hotspots_unchanged": True,
        "mask_has_expected_ten_positions": True,
        "saved_nll_reproduces_scores": True,
        "official_score_function_agrees": True,
        "hard_rule_negative_controls_rejected": True,
        "pareto_orientation_and_eligibility_checked": True,
        "same_random_orders_confirmed_for_batch_comparison": True,
        "batch_size_invariance_max_absolute_error": maximum_error,
        "structural_or_experimental_validation": False}
    (output / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
