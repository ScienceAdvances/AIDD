"""Prepare an isolated copy of the official mouse MDM2-p53 exercise."""
import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/peptide/HADDOCK3-protein-peptide"
AA = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))

def inspect_pdb(path):
    models = {}
    model_id = 1
    for line in path.read_text().splitlines():
        if line.startswith("MODEL "):
            model_id = int(line[10:14])
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
            chain = line[21].strip()
            key = (int(line[22:26]), line[26])
            models.setdefault(model_id, {}).setdefault(chain, {})[key] = line[17:20]
    return [
        {"model": model, "chain": chain, "count": len(residues),
         "first": next(iter(residues))[0], "last": list(residues)[-1][0],
         "sequence": "".join(AA[name] for name in residues.values())}
        for model, chains in models.items() for chain, residues in chains.items()
    ]

def prepare(output):
    source_files = [
        "pdbs/AF_MDM2_26_109.pdb", "pdbs/peptide_ens.pdb", "pdbs/1YCR.pdb",
        "restraints/protein-AR3D-active.act-pass",
        "restraints/peptide-ens-passive.act-pass",
        "restraints/protein-peptide_ambig.tbl",
    ]
    inventories = {rel: inspect_pdb(SOURCE / rel)
                   for rel in source_files if rel.endswith(".pdb")}
    receptor = inventories["pdbs/AF_MDM2_26_109.pdb"]
    peptide = inventories["pdbs/peptide_ens.pdb"]
    if len(receptor) != 1 or any(receptor[0][k] != v for k, v in
            {"chain": "A", "count": 84, "first": 26, "last": 109}.items()):
        raise ValueError("Unexpected receptor chain or residue coverage")
    if len(peptide) != 3 or any(
        row["chain"] != "B" or row["count"] != 15 or row["first"] != 1
        or row["last"] != 15 or row["sequence"] != "SQETFSGLWKLLPPE"
        for row in peptide
    ):
        raise ValueError("Unexpected peptide ensemble")
    active_text = (SOURCE / source_files[3]).read_text().splitlines()
    passive_text = (SOURCE / source_files[4]).read_text().splitlines()
    active = [int(x) for x in active_text[0].split()]
    if active != [54, 55, 58, 59, 62, 67, 72, 73, 93, 94, 100]:
        raise ValueError("Unexpected receptor active residues")
    if active_text[1].strip() or passive_text[0].strip():
        raise ValueError("Unexpected active/passive assignment")
    if [int(x) for x in passive_text[1].split()] != list(range(1, 16)):
        raise ValueError("Unexpected peptide passive residues")
    config = (SOURCE / "workflows/protein_peptide_docking.cfg").read_text()
    output.mkdir(parents=True, exist_ok=False)
    for rel in source_files:
        dest = output / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SOURCE / rel, dest)
    (output / "workflows").mkdir()
    for name, cores, sampling, selected in [
        ("smoke", 2, 12, 12), ("full", 4, 1000, 200)
    ]:
        updated = re.sub(r'^run_dir = .*$', f'run_dir = "runs/{name}"', config, flags=re.M)
        updated = re.sub(r'^ncores = .*$', f'ncores = {cores}', updated, flags=re.M)
        updated = re.sub(r'^sampling = .*$', f'sampling = {sampling}', updated, flags=re.M)
        updated = re.sub(r'^select = .*$', f'select = {selected}', updated, flags=re.M)
        (output / f"workflows/{name}.cfg").write_text(updated)
    hashes = {rel: hashlib.sha256((output / rel).read_bytes()).hexdigest()
              for rel in source_files}
    manifest = {"source": str(SOURCE.relative_to(ROOT)), "inventory": inventories,
                "input_sha256": hashes, "stage": "input_preparation"}
    (output / "preparation.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"output": str(output), "inventory": inventories}, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "workflows/haddock3_practice")
    args = parser.parse_args()
    prepare(args.output.resolve())
