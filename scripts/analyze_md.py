"""Analyze official and locally generated trajectories with explicit selections and references."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
import mdtraj as md
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples/debian13_20261023/md_analysis"
REF = ROOT / "data/md/md-intro-tutorial-main/data/reference"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gmx", default=shutil.which("gmx") or "gmx")
    p.add_argument("--local-run", type=Path, required=True)
    p.add_argument("--output", type=Path, default=OUT)
    args = p.parse_args()
    output = args.output.resolve()
    local_source = args.local_run.resolve()
    output.mkdir(parents=True, exist_ok=True)
    local_label = "local_run"
    summaries = []
    for label, source, topology, trajectory, coords, energy in [
        ("official_1ns", REF, "em_charmm.tpr", "md_charmm.xtc", "em_charmm.gro", "md_charmm.edr"),
        (local_label, local_source, "md.tpr", "md.xtc", "em.gro", "md.edr"),
    ]:
        work = output / label
        work.mkdir(parents=True, exist_ok=True)
        (work / "logs").mkdir(exist_ok=True)
        for filename, target in [(topology, "reference.tpr"), (trajectory, "raw.xtc"), (coords, "reference.gro"), (energy, "md.edr")]:
            dst = work / target
            expected = (source / filename).resolve(strict=True)
            if dst.is_symlink():
                if dst.resolve() != expected:
                    raise ValueError(f"Analysis path points to another input: {dst}")
            elif dst.exists():
                raise FileExistsError(f"Analysis input name already exists: {dst}")
            else:
                dst.symlink_to(expected)
        commands = []

        def run(name, command, stdin=""):
            full = [args.gmx] + command
            commands.append({"name": name, "argv": command, "stdin": stdin})
            with (work / "logs" / (name + ".log")).open("w") as f:
                r = subprocess.run(full, input=stdin, text=True, cwd=work, stdout=f, stderr=subprocess.STDOUT,
                                   env={**os.environ, "GMX_MAXBACKUP": "-1"})
            if r.returncode:
                raise RuntimeError(f"{label}/{name} failed; see its log")

        run("reference", ["trjconv", "-s", "reference.tpr", "-f", "reference.gro", "-o", "reference_protein.pdb", "-pbc", "cluster"], "Protein\nProtein\n")
        run("whole", ["trjconv", "-s", "reference.tpr", "-f", "raw.xtc", "-o", "whole.xtc", "-pbc", "whole"], "Protein\n")
        run("nojump", ["trjconv", "-s", "reference.tpr", "-f", "whole.xtc", "-o", "nojump.xtc", "-pbc", "nojump"], "Protein\n")
        run("center", ["trjconv", "-s", "reference.tpr", "-f", "nojump.xtc", "-o", "protein_centered.xtc", "-center", "-pbc", "mol"], "Protein\nProtein\n")
        run("rmsd_tpr_reference", ["rms", "-s", "reference.tpr", "-f", "protein_centered.xtc", "-o", "rmsd_ca_to_tpr_reference.xvg", "-tu", "ns", "-xvg", "none"], "C-alpha\nC-alpha\n")
        run("rmsd", ["rms", "-s", "reference_protein.pdb", "-f", "protein_centered.xtc", "-o", "rmsd_ca_to_em.xvg", "-tu", "ns", "-xvg", "none"], "C-alpha\nC-alpha\n")
        run("rmsf", ["rmsf", "-s", "reference.tpr", "-f", "protein_centered.xtc", "-o", "rmsf_ca.xvg", "-res", "-xvg", "none"], "C-alpha\n")
        run("energy", ["energy", "-f", "md.edr", "-o", "thermodynamics.xvg", "-xvg", "none"], "Temperature\nPressure\nDensity\n0\n")
        (work / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        # XTC readers may reject non-ASCII paths; pass relative ASCII names from the data folder.
        before = Path.cwd()
        os.chdir(work)
        try:
            t = md.load("protein_centered.xtc", top="reference_protein.pdb")
            reference = md.load("reference_protein.pdb")
        finally:
            os.chdir(before)
        ca = t.top.select("name CA")
        masses = np.array([atom.element.mass for atom in t.top.atoms])
        rg = md.compute_rg(t, masses=masses)
        rmsd_first = md.rmsd(t, t, 0, atom_indices=ca, parallel=False)
        gmx_rmsd = np.loadtxt(work / "rmsd_ca_to_em.xvg", ndmin=2)
        assert len(gmx_rmsd) == t.n_frames
        assert np.isfinite(rg).all() and np.isfinite(rmsd_first).all()
        assert np.max(np.abs(gmx_rmsd[:, 0] - t.time / 1000)) < 1e-5
        independent_rmsd = md.rmsd(t, reference, atom_indices=ca, parallel=False)
        rmsd_error = float(np.max(np.abs(gmx_rmsd[:, 1] - independent_rmsd)))
        assert rmsd_error < 2e-5
        series = pd.DataFrame({"time_ps": t.time, "time_ns": t.time / 1000,
                               "ca_rmsd_to_minimized_reference_nm": gmx_rmsd[:, 1],
                               "ca_rmsd_to_first_md_frame_nm": rmsd_first, "protein_mass_weighted_rg_nm": rg})
        series.to_csv(work / "structural_metrics.csv", index=False)
        therm = np.loadtxt(work / "thermodynamics.xvg", ndmin=2)
        assert therm.shape[1] == 4
        pd.DataFrame(therm, columns=["time_ps", "temperature_K", "pressure_bar", "density_kg_m3"]).to_csv(work / "thermodynamics.csv", index=False)
        s = {"source": label, "input_directory": str(source), "frames": t.n_frames, "protein_atoms": t.n_atoms, "ca_atoms": len(ca),
             "start_ps": float(t.time[0]), "end_ps": float(t.time[-1]),
             "rmsd_to_em_mean_nm": float(gmx_rmsd[:, 1].mean()), "rmsd_to_first_mean_nm": float(rmsd_first.mean()),
             "rg_mean_nm": float(rg.mean()), "temperature_mean_K": float(therm[:, 1].mean()),
             "pressure_mean_bar": float(therm[:, 2].mean()), "density_mean_kg_m3": float(therm[:, 3].mean()),
             "independent_rmsd_max_error_nm": rmsd_error,
             "relevant_sampling_convergence_claim": False}
        summaries.append(s)
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7), constrained_layout=True)
        axes[0, 0].plot(series.time_ns, series.ca_rmsd_to_minimized_reference_nm, label="Reference: EM structure")
        axes[0, 0].plot(series.time_ns, series.ca_rmsd_to_first_md_frame_nm, label="Reference: first MD frame")
        axes[0, 0].set(xlabel="Time (ns)", ylabel="CA RMSD (nm)"); axes[0, 0].legend(fontsize=8)
        axes[0, 1].plot(series.time_ns, rg, color="#a96036")
        axes[0, 1].set(xlabel="Time (ns)", ylabel="Protein mass-weighted Rg (nm)")
        axes[1, 0].plot(therm[:, 0] / 1000, therm[:, 1], color="#19627e")
        axes[1, 0].axhline(300, color="#7c7c7c", ls="--")
        axes[1, 0].set(xlabel="Time (ns)", ylabel="Temperature (K)")
        axes[1, 1].plot(therm[:, 0] / 1000, therm[:, 3], color="#50806b")
        axes[1, 1].set(xlabel="Time (ns)", ylabel="Density (kg/m3)")
        for ax in axes.ravel():ax.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"Factor Xa in water: {label}")
        fig.savefig(work / "trajectory_analysis.png", dpi=180); plt.close(fig)
        print(json.dumps(s))
    (output / "analysis_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")


if __name__ == "__main__":
    main()
