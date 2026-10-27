"""Run official ProteinMPNN on the experimental 1YCR complex with audited masks."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT / "data/aidd/ProteinMPNN"
COMMIT = "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", default="runs/02_proteinmpnn")
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = p.parse_args()
    import torch
    device = "cuda" if args.device == "auto" and torch.cuda.is_available() else args.device
    if device == "auto":
        device = "cpu"
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    out = (ROOT / args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "pdbs").mkdir()
    (out / "logs").mkdir()
    pdb = ROOT / "data/peptide/1YCR_original.pdb"
    shutil.copy2(pdb, out / "pdbs/1YCR.pdb")
    records = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
            key = (line[21], int(line[22:26]), line[26])
            records.setdefault(key, line[17:20].strip())
    b = [(key, aa) for key, aa in records.items() if key[0] == "B"]
    fixed = [i + 1 for i, (key, aa) in enumerate(b) if key[1] in (19, 23, 26)]
    assert fixed == [3, 7, 10]
    assert [b[i - 1][1] for i in fixed] == ["PHE", "TRP", "LEU"]
    with (out / "numbering.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["parsed_position", "pdb_chain", "pdb_residue", "residue", "fixed"])
        for i, (key, aa) in enumerate(b):
            w.writerow([i + 1, key[0], key[1], aa, i + 1 in fixed])
    commands = []
    started = time.monotonic()

    def run(name, script, argv):
        command = [sys.executable, str(REPO / script)] + argv
        commands.append({"step": name, "argv": command})
        with (out / "logs" / (name + ".log")).open("w") as f:
            r = subprocess.run(command, cwd=REPO, stdout=f, stderr=subprocess.STDOUT,
                               env={**os.environ, "OMP_NUM_THREADS": "4", "MKL_NUM_THREADS": "4",
                                    **({"CUDA_VISIBLE_DEVICES": ""} if device == "cpu" else {})})
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        if r.returncode:
            raise RuntimeError(f"{name} failed; inspect its log")

    parsed = str(out / "parsed.jsonl")
    chains = str(out / "chains.jsonl")
    positions = str(out / "fixed_positions.jsonl")
    run("parse", "helper_scripts/parse_multiple_chains.py", ["--input_path", str(out / "pdbs"), "--output_path", parsed])
    parsed_data = json.loads(Path(parsed).read_text().splitlines()[0])
    assert parsed_data["seq_chain_B"] == "ETFSDLWKLLPEN"
    assert len(parsed_data["seq_chain_A"]) == 85
    run("chains", "helper_scripts/assign_fixed_chains.py", ["--input_path", parsed, "--output_path", chains, "--chain_list", "B"])
    run("fixed", "helper_scripts/make_fixed_positions_dict.py", ["--input_path", parsed, "--output_path", positions, "--chain_list", "B", "--position_list", "3 7 10"])
    run("design", "protein_mpnn_run.py", ["--jsonl_path", parsed, "--chain_id_jsonl", chains,
        "--fixed_positions_jsonl", positions, "--out_folder", str(out / "designs"),
        "--model_name", "v_48_020", "--num_seq_per_target", "2", "--batch_size", "1",
        "--sampling_temp", "0.1 0.3 0.5", "--seed", "20261010", "--save_score", "1", "--save_probs", "1"])
    lines = (out / "designs/seqs/1YCR.fa").read_text().splitlines()
    native = lines[1]
    candidates = []
    for i in range(2, len(lines), 2):
        header, seq = lines[i], lines[i + 1]
        assert len(seq) == 13 and all(seq[j - 1] == native[j - 1] for j in fixed)
        fields = dict(re.findall(r"(T|sample|score|global_score|seq_recovery)=([\d.\-]+)", header))
        mutable = [j for j in range(13) if j + 1 not in fixed]
        identity = sum(seq[j] == native[j] for j in mutable) / len(mutable)
        candidates.append({"candidate": f"T{fields['T']}_sample{fields['sample']}", "sequence": seq,
                           "temperature": float(fields["T"]), "score": float(fields["score"]),
                           "global_score": float(fields["global_score"]),
                           "mutable_identity": identity, "official_seq_recovery": float(fields["seq_recovery"])})
    with (out / "candidate_summary.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(candidates[0]))
        w.writeheader(); w.writerows(candidates)
    summary = {"completed": True, "upstream_commit": COMMIT, "device": device,
               "executed_at": __import__("datetime").datetime.now().astimezone().isoformat(),
               "pdb": "1YCR", "native_chain_B": native, "fixed_chain": "A",
               "fixed_positions_B": fixed, "mutable_positions": 10,
               "model": "v_48_020", "seed": 20261010, "candidates": candidates,
               "elapsed_seconds": time.monotonic() - started,
               "torch_version": __import__("torch").__version__,
               "checkpoint_sha256": hashlib.sha256((REPO / "vanilla_model_weights/v_48_020.pt").read_bytes()).hexdigest(),
               "experimental_validation": False}
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
