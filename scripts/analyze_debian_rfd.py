"""Summarize native RFD2/Chai metrics using the actual compiled column names."""
from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'examples/debian13_20261023'
SOURCE = ROOT / 'examples/debian13_20261023/rfdiffusion2_outputs/compiled_metrics.csv'


def main():
    frame = pd.read_csv(SOURCE)
    assert len(frame) == 8
    rows=[]
    for index, row in frame.iterrows():
        model=int(row['pocket_aligned_total_ligand_rmsd_argmin'])
        values=dict(candidate_row=index+1, source_csv_line=index+2,
            model_with_minimum_ligand_rmsd=model,
            motif_rmsd_to_input_A=float(row['motif_rmsd_to_input']),
            idealized_motif_constellation_rmsd_A=float(row['metrics.IdealizedResidueRMSD.rmsd_constellation']),
            minimum_pocket_aligned_ligand_rmsd_A=float(row['pocket_aligned_total_ligand_rmsd_min']),
            complex_plddt_for_that_model=float(row[f'complex_plddt_chaimodel_{model}_sym_resolved']),
            backbone_aligned_motif_atom_rmsd_to_reference_A=float(row[f'backbone_aligned_contigatom_rmsd_chai_ref_chaimodel_{model}_sym_resolved']))
        rows.append(values)
    result=pd.DataFrame(rows)
    result.to_csv(OUT/'rfdiffusion2_candidate_summary.csv',index=False)
    summary=dict(completed=True, backbones=1, sequence_candidates=8, chai_predictions=40,
        source=str(SOURCE), compiled_shape=list(frame.shape), optional_esm_embeddings=False,
        minimum_ligand_rmsd_range_A=[float(result.minimum_pocket_aligned_ligand_rmsd_A.min()),float(result.minimum_pocket_aligned_ligand_rmsd_A.max())],
        rows=rows, experimental_activity_evaluated=False)
    (OUT/'rfdiffusion2_analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(7,4.5),constrained_layout=True)
    ax.scatter(result.complex_plddt_for_that_model,result.minimum_pocket_aligned_ligand_rmsd_A,s=60,color='#327a88')
    for _,r in result.iterrows():
        ax.annotate(str(int(r.candidate_row)),(r.complex_plddt_for_that_model,r.minimum_pocket_aligned_ligand_rmsd_A),xytext=(5,4),textcoords='offset points')
    ax.set(xlabel='Complex pLDDT (0–1) for that sample',ylabel='Minimum ligand RMSD over 5 samples (A)',
        title='Native Chai refolding: confidence and ligand geometry')
    fig.savefig(OUT/'rfdiffusion2_confidence_and_ligand_geometry.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))



if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    SOURCE=args.source.resolve()
    OUT=args.output.resolve()
    OUT.mkdir(parents=True, exist_ok=False)
    main()
