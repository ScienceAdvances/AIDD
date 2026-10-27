"""Reproduce the official Factor Xa setup; default to a labeled short CPU run."""
from pathlib import Path
import argparse
import json
import os
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
OFFICIAL = ROOT / "data/md/md-intro-tutorial-main/data"


def set_parameter(text, key, value):
    pattern = rf"(?mi)^\s*{re.escape(key)}\s*=.*$"
    line = f"{key} = {value}"
    return re.sub(pattern, line, text) if re.search(pattern, text) else text + "\n" + line + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gmx", default=shutil.which("gmx") or "gmx")
    p.add_argument("--full", action="store_true")
    p.add_argument("--output", default="runs/03_md")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--gpu", action="store_true", help="Require CUDA nonbonded offload for mdrun")
    args = p.parse_args()
    output = ROOT / args.output
    output.mkdir(parents=True, exist_ok=False)
    (output / "logs").mkdir()
    commands = []
    started = time.monotonic()

    def run(name, argv, stdin=""):
        if args.gpu and argv[0] == "mdrun":
            argv = argv + ["-nb", "gpu", "-gpu_id", "0"]
        full = [args.gmx] + argv
        commands.append({"step": name, "argv": full, "stdin": stdin})
        print(name, flush=True)
        with (output / "logs" / f"{name}.log").open("w") as f:
            result = subprocess.run(full, input=stdin, text=True, cwd=output, stdout=f, stderr=subprocess.STDOUT,
                                    env={**os.environ, "GMX_MAXBACKUP": "-1"})
        if result.returncode:
            raise RuntimeError(f"{name} failed; inspect {output / 'logs' / (name + '.log')}")
        (output / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")

    pdb = (OFFICIAL / "input/1fjs.pdb").read_text()
    cleaned = "\n".join(line for line in pdb.splitlines() if not line.startswith(("HETATM", "CONECT"))) + "\n"
    (output / "protein.pdb").write_text(cleaned)
    run("version", ["--version"])
    run("topology", ["pdb2gmx", "-f", "protein.pdb", "-o", "processed.gro", "-p", "topol.top", "-ff", "charmm27", "-water", "tip3p", "-ignh"], "\n" * 16)
    run("box", ["editconf", "-f", "processed.gro", "-o", "box.gro", "-c", "-d", "1.0", "-bt", "dodecahedron"])
    run("solvate", ["solvate", "-cp", "box.gro", "-cs", "spc216.gro", "-o", "solv.gro", "-p", "topol.top"])
    ions = (OFFICIAL / "input/emin-charmm.mdp").read_text()
    ions = set_parameter(ions, "coulombtype", "Cut-off")
    ions = set_parameter(ions, "nsteps", 0)
    (output / "ions.mdp").write_text(ions)
    run("ions_input", ["grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "topol.top", "-o", "ions.tpr"])
    run("ions", ["genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top", "-pname", "NA", "-nname", "CL", "-neutral", "-conc", "0.15", "-seed", "20261010"], "SOL\n")
    durations = {"nvt": 10, "npt": 10, "md": 20}
    for stage, filename in [("em", "emin-charmm.mdp"), ("nvt", "nvt-charmm.mdp"), ("npt", "npt-charmm.mdp"), ("md", "md-charmm.mdp")]:
        params = (OFFICIAL / "input" / filename).read_text()
        if not args.full and stage != "em":
            dt = float(re.search(r"(?mi)^\s*dt\s*=\s*([\d.]+)", params).group(1))
            params = set_parameter(params, "nsteps", round(durations[stage] / dt))
            for key in ["nstlog", "nstenergy", "nstxout-compressed"]:
                params = set_parameter(params, key, 50)
        if stage == "nvt":
            params = set_parameter(params, "gen_seed", 20261010)
        (output / f"{stage}.mdp").write_text(params)
    run("em_input", ["grompp", "-f", "em.mdp", "-c", "solv_ions.gro", "-p", "topol.top", "-o", "em.tpr"])
    run("em", ["mdrun", "-deffnm", "em", "-ntmpi", "1", "-ntomp", str(args.threads)])
    run("nvt_input", ["grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro", "-p", "topol.top", "-o", "nvt.tpr"])
    run("nvt", ["mdrun", "-deffnm", "nvt", "-ntmpi", "1", "-ntomp", str(args.threads)])
    run("npt_input", ["grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro", "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr"])
    run("npt", ["mdrun", "-deffnm", "npt", "-ntmpi", "1", "-ntomp", str(args.threads)])
    run("md_input", ["grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt", "-p", "topol.top", "-o", "md.tpr"])
    run("md", ["mdrun", "-deffnm", "md", "-ntmpi", "1", "-ntomp", str(args.threads)])
    gro = (output / "md.gro").read_text().splitlines()
    summary = {"complete": True, "mode": "official_durations" if args.full else "short_workflow_demonstration",
               "stage_ps": "see saved mdp files" if args.full else durations, "atoms": int(gro[1]),
               "force_field": "charmm27", "water": "tip3p", "ion_concentration_M": 0.15,
               "gpu_nonbonded_required": args.gpu,
               "elapsed_seconds": time.monotonic() - started, "command_count": len(commands)}
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
