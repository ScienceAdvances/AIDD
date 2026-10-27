"""Recalculate CA RMSD from the distributed protein-only XTC using MDTraj."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import mdtraj as md
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    source=ROOT/'examples/debian13_20261023/md_analysis/local_run'
    trajectory=md.load(str(source/'protein_centered.xtc'),top=str(source/'reference_protein.pdb'))
    reference=md.load(str(source/'reference_protein.pdb'))
    ca=trajectory.top.select('name CA');rmsd=md.rmsd(trajectory,reference,atom_indices=ca,parallel=False)
    with (source/'structural_metrics.csv').open() as f:rows=list(csv.DictReader(f))
    expected=np.array([float(r['ca_rmsd_to_minimized_reference_nm']) for r in rows]);error=float(np.max(np.abs(rmsd-expected)))
    if trajectory.n_frames!=201 or error>1e-5:raise ValueError(f'Trajectory/CSV mismatch: {error}')
    a.output.mkdir(parents=True,exist_ok=False)
    result=dict(passed=True,frames=201,ca_atoms=len(ca),max_RMSD_error_nm=error,mean_RMSD_nm=float(rmsd.mean()),trajectory_contains_water=False)
    (a.output/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
