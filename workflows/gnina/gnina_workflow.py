"""Auditable GNINA preparation, execution, score extraction and pose validation."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request


def download(url, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    request = urllib.request.Request(url, headers={"User-Agent": "GNINA-Tutorial/1.0"})
    with urllib.request.urlopen(request, timeout=90) as response:
        content = response.read()
    if not content:
        raise ValueError(f"Empty response: {url}")
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_bytes(content)
    temporary.replace(path)
    return {"url": url, "file": str(path), "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest()}


def chemistry():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import rdkit
    except ImportError as exc:
        raise RuntimeError("Install rdkit in the selected Python environment") from exc
    return Chem, AllChem, rdkit.__version__


def read_sdf(path):
    Chem, _, _ = chemistry()
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if not molecules or any(m is None for m in molecules):
        raise ValueError(f"Empty or invalid SDF record: {path}")
    for mol in molecules:
        if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
            raise ValueError(f"Expected one 3D conformer per SDF record: {path}")
        xyz = mol.GetConformer().GetPositions()
        if not all(math.isfinite(float(v)) for row in xyz for v in row):
            raise ValueError(f"Non-finite coordinates: {path}")
    return molecules


def example(args):
    Chem, AllChem, version = chemistry()
    out = args.work.resolve()
    out.mkdir(parents=True, exist_ok=False)
    base = "https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/basic_docking/data/"
    sources = [download(base + "1iep_receptorH.pdb", out / "receptor.pdb"),
               download(base + "1iep_ligand.sdf", out / "reference.sdf")]
    native = read_sdf(out / "reference.sdf")[0]
    # Generate an unbound starting geometry while preserving source chemistry.
    ligand = Chem.AddHs(Chem.RemoveHs(native))
    ligand.RemoveAllConformers()
    settings = AllChem.ETKDGv3()
    settings.randomSeed = 20261011
    if AllChem.EmbedMolecule(ligand, settings) != 0:
        raise RuntimeError("ETKDG embedding failed")
    if not AllChem.MMFFHasAllMoleculeParams(ligand):
        raise RuntimeError("MMFF parameters are unavailable")
    optimization = AllChem.MMFFOptimizeMolecule(ligand, maxIters=2000)
    if optimization != 0:
        raise RuntimeError(f"MMFF optimization did not converge: {optimization}")
    ligand.SetProp("_Name", "imatinib_official_source_state")
    writer = Chem.SDWriter(str(out / "ligand_input.sdf"))
    writer.write(ligand); writer.close()
    coordinates = native.GetConformer().GetPositions()
    lower, upper = coordinates.min(axis=0), coordinates.max(axis=0)
    center = (lower + upper) / 2
    size = upper - lower + 8.0
    native_heavy = Chem.RemoveHs(native)
    record = {"sources": sources, "rdkit_version": version,
              "heavy_atoms": native_heavy.GetNumAtoms(),
              "source_formal_charge": Chem.GetFormalCharge(native),
              "source_isomeric_smiles": Chem.MolToSmiles(native_heavy),
              "generated_conformer_seed": 20261011,
              "mmff_converged": optimization == 0,
              "reference_center_A": center.tolist(),
              "reference_extent_plus_8A": size.tolist(),
              "box_note": "Reference bounds only; GNINA autobox_extend may enlarge the actual box",
              "ph_state_enumeration_completed": False,
              "gnina_docking_completed": False}
    (out / "input_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


def pubchem(args):
    if args.cid <= 0:
        raise ValueError("CID must be positive")
    record = download(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{args.cid}/SDF?record_type=3d",
        args.output)
    read_sdf(args.output)
    Path(str(args.output) + ".source.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record))


def check(args):
    Chem, _, _ = chemistry()
    rows = []
    for path in args.sdf:
        for index, molecule in enumerate(read_sdf(path), 1):
            rows.append({"file": str(path), "record": index,
                         "heavy_atoms": molecule.GetNumHeavyAtoms(),
                         "formal_charge": Chem.GetFormalCharge(molecule),
                         "smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule))})
    print(json.dumps(rows, indent=2))


def jobs_from_file(path):
    names = ["job_id", "receptor", "ligand", "reference", "seed"]
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != names:
            raise ValueError(f"TSV header must be: {names}")
        rows = list(reader)
    if not rows:
        raise ValueError("No jobs supplied")
    ids = set()
    for row in rows:
        key = row["job_id"]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", key) or key in ids:
            raise ValueError(f"Unsafe or duplicate job_id: {key}")
        ids.add(key)
        row["seed"] = int(row["seed"])
        if row["seed"] <= 0:
            raise ValueError("Use an explicit positive seed")
        for name in ["receptor", "ligand", "reference"]:
            target = Path(row[name])
            row[name] = str((target if target.is_absolute() else path.parent / target).resolve())
            if not Path(row[name]).is_file():
                raise FileNotFoundError(row[name])
        if row["receptor"] == row["reference"]:
            raise ValueError("Pocket workflow requires a ligand reference, not the receptor")
        if len(read_sdf(row["reference"])) != 1:
            raise ValueError("Use one selected reference ligand for one pocket")
        read_sdf(row["ligand"])
    return rows


def run(args):
    if args.cpu < 1 or args.exhaustiveness < 1 or args.num_modes < 1:
        raise ValueError("CPU count, exhaustiveness and num_modes must be positive")
    rows = jobs_from_file(args.jobs.resolve())
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite run directory: {out}")
    binary = shutil.which(args.gnina) if not Path(args.gnina).is_file() else str(Path(args.gnina).resolve())
    if not binary and not args.dry_run:
        raise FileNotFoundError("GNINA is unavailable; use --dry-run to inspect commands")
    command_binary = binary or args.gnina
    out.mkdir(parents=True)
    version = "not_executed_dry_run"
    if not args.dry_run:
        version = subprocess.check_output([command_binary, "--version"], text=True, stderr=subprocess.STDOUT).strip()
        (out / "gnina_version.txt").write_text(version + "\n")
    records = []
    failures = 0
    # Execute sequentially by design: one job uses the declared CPU/GPU budget.
    for row in rows:
        key = row["job_id"]
        pose = out / f"{key}.sdf"
        command = [command_binary, "--receptor", row["receptor"], "--ligand", row["ligand"],
                   "--autobox_ligand", row["reference"], "--autobox_add", "4",
                   "--out", str(pose), "--cpu", str(args.cpu),
                   "--scoring", "vina", "--cnn_scoring", args.cnn_scoring,
                   "--pose_sort_order", "Energy" if args.cnn_scoring == "none" else "CNNscore",
                   "--seed", str(row["seed"]), "--exhaustiveness", str(args.exhaustiveness),
                   "--num_modes", str(args.num_modes)]
        command += ["--no_gpu"] if args.device is None else ["--device", str(args.device)]
        record = {"job": row, "argv": command, "version": version,
                  "input_sha256": {name: hashlib.sha256(Path(row[name]).read_bytes()).hexdigest()
                                   for name in ["receptor", "ligand", "reference"]},
                  "dry_run": args.dry_run}
        if args.dry_run:
            print(json.dumps(command))
        else:
            start = time.monotonic()
            with (out / f"{key}.log").open("w") as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                        env={**os.environ, "OMP_NUM_THREADS": str(args.cpu),
                                             "MKL_NUM_THREADS": str(args.cpu)})
            record["returncode"] = result.returncode
            record["elapsed_seconds"] = time.monotonic() - start
            record["success"] = False
            if result.returncode == 0:
                try:
                    record["output_records"] = len(read_sdf(pose))
                    record["success"] = True
                except (ValueError, FileNotFoundError) as exc:
                    record["error"] = str(exc)
            failures += not record["success"]
        records.append(record)
        (out / "run_manifest.json").write_text(json.dumps(records, indent=2) + "\n")
    if failures:
        raise RuntimeError(f"{failures} jobs failed; inspect run_manifest.json and logs")


def property_value(molecule, aliases):
    for name in aliases:
        if molecule.HasProp(name):
            raw = molecule.GetProp(name).strip()
            try:
                number = float(raw)
            except ValueError as exc:
                raise ValueError(f"Non-numeric score property {name}: {raw}") from exc
            if not math.isfinite(number):
                raise ValueError(f"Non-finite score property {name}")
            return number
    return None


def summarize(args):
    Chem, _, _ = chemistry()
    run_path = args.run.resolve()
    manifest = json.loads((run_path / "run_manifest.json").read_text())
    if any(item["dry_run"] for item in manifest):
        raise ValueError("A dry run has no docking results to summarize")
    all_rows = []
    failures = []
    for record in manifest:
        if not record["success"]:
            failures.append(record["job"]["job_id"])
            continue
        job = record["job"]
        for rank, molecule in enumerate(read_sdf(run_path / (job["job_id"] + ".sdf")), 1):
            all_rows.append({"job_id": job["job_id"], "seed": job["seed"],
                             "receptor": job["receptor"], "ligand": job["ligand"],
                             "output_rank": rank,
                             "molecule_name": molecule.GetProp("_Name") if molecule.HasProp("_Name") else "",
                             "vina_score_kcal_mol": property_value(molecule, ["minimizedAffinity", "Affinity"]),
                             "CNNscore": property_value(molecule, ["CNNscore"]),
                             "CNNaffinity_pK": property_value(molecule, ["CNNaffinity"]),
                             "CNN_VS": property_value(molecule, ["CNN_VS"]),
                             "heavy_atoms": molecule.GetNumHeavyAtoms(),
                             "source_sdf": str(run_path / (job["job_id"] + ".sdf"))})
    if not all_rows:
        raise ValueError("No successful output molecules")
    if args.output.exists():
        raise FileExistsError(args.output)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0]))
        writer.writeheader(); writer.writerows(all_rows)
    print(json.dumps({"pose_rows": len(all_rows), "failed_jobs": failures,
                      "missing_properties_preserved_as_empty": True}))


def rmsd(args):
    Chem, AllChem, _ = chemistry()
    references = read_sdf(args.reference)
    if len(references) != 1:
        raise ValueError("RMSD requires exactly one reference ligand")
    reference = Chem.RemoveHs(references[0])
    rows = []
    identity = Chem.MolToSmiles(reference)
    for rank, molecule in enumerate(read_sdf(args.poses), 1):
        probe = Chem.RemoveHs(molecule)
        if Chem.MolToSmiles(probe) != identity:
            raise ValueError(f"Pose {rank} differs in chemistry; resolve mapping before RMSD")
        # CalcRMS keeps coordinates in the receptor frame; do not fit the ligand itself.
        value = AllChem.CalcRMS(probe, reference, maxMatches=100000)
        rows.append({"output_rank": rank, "heavy_atom_rmsd_A": value,
                     "within_2A": value <= 2.0})
    if args.output.exists():
        raise FileExistsError(args.output)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"top1_rmsd_A": rows[0]["heavy_atom_rmsd_A"],
                      "best_sampled_rmsd_A": min(r["heavy_atom_rmsd_A"] for r in rows),
                      "ligand_fitting_applied": False}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    a = sub.add_parser("fetch-example"); a.add_argument("--work", type=Path, required=True); a.set_defaults(func=example)
    a = sub.add_parser("fetch-pubchem"); a.add_argument("--cid", type=int, required=True); a.add_argument("--output", type=Path, required=True); a.set_defaults(func=pubchem)
    a = sub.add_parser("check"); a.add_argument("--sdf", type=Path, nargs="+", required=True); a.set_defaults(func=check)
    a = sub.add_parser("run"); a.add_argument("--jobs", type=Path, required=True); a.add_argument("--output", type=Path, required=True)
    a.add_argument("--gnina", default="gnina"); a.add_argument("--cpu", type=int, default=4)
    a.add_argument("--device", type=int); a.add_argument("--exhaustiveness", type=int, default=32)
    a.add_argument("--num-modes", type=int, default=20)
    a.add_argument("--cnn-scoring", choices=["none", "rescore"], default="rescore")
    a.add_argument("--dry-run", action="store_true"); a.set_defaults(func=run)
    a = sub.add_parser("summarize"); a.add_argument("--run", type=Path, required=True); a.add_argument("--output", type=Path, required=True); a.set_defaults(func=summarize)
    a = sub.add_parser("rmsd"); a.add_argument("--reference", type=Path, required=True); a.add_argument("--poses", type=Path, required=True); a.add_argument("--output", type=Path, required=True); a.set_defaults(func=rmsd)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
