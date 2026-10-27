"""Recalculate lightweight Debian evidence without executing a design model."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
EX=ROOT/'examples/debian13_20261023'

def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))

def audit_sequence(output):
    directory=EX/'sequence_debian_only';candidates=rows(directory/'candidates.csv')
    with np.load(directory/'score_audit.npz',allow_pickle=False) as saved:
        mask=saved['score_mask'];scores=saved['scores'];nll=saved['per_position_nll']
        recalculated=(nll*mask).sum(axis=2)/mask.sum()
        if not np.allclose(recalculated,scores,atol=1e-6):raise ValueError('NLL arithmetic mismatch')
        maximum=float(np.max(np.abs(recalculated-scores)))
        if int(mask.sum())!=10:raise ValueError('Wrong scoring mask')
        if not np.allclose(scores.mean(axis=1),[float(r['mpnn_nll_mean']) for r in candidates],atol=1e-6):raise ValueError('CSV mismatch')
    native='ETFSDLWKLLPEN'
    for r in candidates:
        seq=r['sequence']
        valid=(len(seq)==13 and all(seq[i-1]==native[i-1] for i in [3,7,10]) and sum(a!=b for a,b in zip(seq,native))<=2 and 'C' not in seq and not any(seq[i-1]=='P' and native[i-1]!='P' for i in range(4,10)))
        if valid!=(r['hard_pass']=='True'):raise ValueError('Hard-rule mismatch')
    fields=['mpnn_nll_mean','gravy','mutation_count']
    eligible=[r for r in candidates if r['hard_pass']=='True']
    pareto=[r for r in eligible if not any(all(float(o[k])<=float(r[k]) for k in fields) and any(float(o[k])<float(r[k]) for k in fields) for o in eligible)]
    if {r['candidate_id'] for r in pareto}!={r['candidate_id'] for r in candidates if r['pareto']=='True'}:raise ValueError('Pareto mismatch')
    selected=[r for r in candidates if r['shortlisted']=='True']
    if (len(candidates),len(eligible),len(pareto),len(selected))!=(197,176,5,4):raise ValueError('Unexpected funnel')
    result=dict(passed=True,candidates=197,eligible=176,pareto=5,shortlist=4,score_mask_positions=10,max_arithmetic_error=maximum,model_inference_executed=False)
    (output/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    with (output/'pareto.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(candidates[0]));w.writeheader();w.writerows(pareto)
    return result

def analyze_md(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    source=EX/'md_analysis/local_run';data=rows(source/'structural_metrics.csv')
    # Preserve original column definitions; identify only the actual headers.
    keys=list(data[0]);print('Structural columns:',keys)
    numeric={k:np.array([float(r[k]) for r in data]) for k in keys}
    rmsd_key='ca_rmsd_to_minimized_reference_nm'
    time_key=next(k for k in keys if 'time' in k.lower())
    reference=next(r for r in json.loads((EX/'md_analysis/analysis_summary.json').read_text()) if r['source']=='local_run')
    mean=float(numeric[rmsd_key].mean())
    if len(data)!=201 or abs(mean-reference['rmsd_to_em_mean_nm'])>1e-6:raise ValueError('MD summary mismatch')
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained');ax.plot(numeric[time_key],numeric[rmsd_key]);ax.set(xlabel=time_key,ylabel=rmsd_key,title='Debian Factor Xa: saved 1 ns trajectory');fig.savefig(output/'rmsd.png',dpi=160);plt.close(fig)
    result=dict(passed=True,frames=len(data),mean_rmsd_to_em_nm=mean,new_dynamics_executed=False,time_column=time_key)
    (output/'verification.json').write_text(json.dumps(result,indent=2)+'\n');return result

def pyrosetta(output):
    data=rows(EX/'pyrosetta/scores.csv')
    if len(data)!=4:raise ValueError('Expected WT plus three variants')
    # Check that side-chain-only refinement preserves CA coordinates.
    from Bio.PDB import PDBParser
    parser=PDBParser(QUIET=True)
    ref=parser.get_structure('ref',ROOT/'data/peptide/1YCR_original.pdb')[0]
    ca={(c.id,r.id[1]):r['CA'].coord for c in ref for r in c if r.id[0]==' ' and 'CA' in r}
    errors={}
    for r in data:
        model=parser.get_structure(r['candidate'],EX/'pyrosetta'/f"{r['candidate']}.pdb")[0]
        diffs=[np.linalg.norm(atom['CA'].coord-ca[(c.id,atom.id[1])]) for c in model for atom in c if atom.id[0]==' ' and 'CA' in atom]
        errors[r['candidate']]=float(max(diffs))
    if max(errors.values())>0.005:raise ValueError('Backbone changed')
    result=dict(passed=True,candidates=4,max_CA_displacement_A=errors,new_refinement_executed=False,energy_unit='REU',experimental_affinity_measured=False)
    (output/'verification.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('chapter',choices=['03','06','10']);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    print(json.dumps({'03':analyze_md,'06':audit_sequence,'10':pyrosetta}[a.chapter](a.output),indent=2))
