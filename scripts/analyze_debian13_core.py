"""Audit fresh HADDOCK and GNINA outputs without reusing historical scores."""
from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'examples/debian13_20261023'


SOURCE = ROOT / 'examples/debian13_20261023'

def main():
    dest = OUT / 'haddock_analysis'
    dest.mkdir(exist_ok=True)
    rows = []
    for run in ['smoke', 'full']:
        for stage in ['02_caprieval', '04_caprieval', '06_caprieval', '13_caprieval']:
            table_path = ROOT / 'examples/haddock_stages' / run / stage / 'capri_ss.tsv'
            if not table_path.exists():
                continue
            frame = pd.read_csv(table_path, sep='\t', comment='#').sort_values('score').reset_index(drop=True)
            best = frame.loc[frame.dockq.idxmax()]
            recalculated = (frame.fnat + 1/(1+(frame.irmsd/1.5)**2) + 1/(1+(frame.lrmsd/8.5)**2))/3
            error = float((recalculated - frame.dockq).abs().max())
            assert error < 0.001
            top = frame.iloc[0]
            rows.append(dict(run=run, stage=stage, models=len(frame), top_score=float(top.score),
                top_score_dockq=float(top.dockq), best_dockq=float(best.dockq),
                best_dockq_score_rank=int(frame.dockq.idxmax())+1, dockq_rounding_error=error,
                table=str(table_path)))
            if stage == '13_caprieval':
                score_error = float((frame.vdw + 0.2*frame.elec + frame.desolv + 0.1*frame.air - frame.score).abs().max())
                assert score_error < 0.002
                frame.to_csv(dest / f'{run}_final_models.csv', index=False)
    pd.DataFrame(rows).to_csv(dest / 'stage_summary.csv', index=False)
    (dest / 'validation.json').write_text(json.dumps(dict(passed=True, rows=rows,
        reference='human 1YCR used to evaluate mouse docking inputs; homology caveat applies'), indent=2)+'\n')
    fig, ax = plt.subplots(figsize=(7.5,4.2), constrained_layout=True)
    for run in ['smoke','full']:
        sub = pd.DataFrame([r for r in rows if r['run']==run])
        ax.plot(sub.stage, sub.best_dockq, 'o-', label=f'{run}: best sampled')
        ax.plot(sub.stage, sub.top_score_dockq, 's--', label=f'{run}: top score')
    ax.set(ylabel='DockQ vs homologous 1YCR', xlabel='Evaluation stage', title='Fresh Debian HADDOCK runs')
    ax.legend(fontsize=8); fig.savefig(dest/'sampling_and_ranking.png',dpi=180);plt.close(fig)
    gnina=[]
    for mode in ['rescore','none']:
        directory=SOURCE/f'gnina_{mode}_gpu'
        if not (directory/'scores.csv').exists(): continue
        scores=pd.read_csv(directory/'scores.csv')
        for f in sorted(directory.glob('*_rmsd.csv')):
            rmsd=pd.read_csv(f);job=f.name.removesuffix('_rmsd.csv')
            top=scores[(scores.job_id==job)&(scores.output_rank==1)].iloc[0]
            gnina.append(dict(mode=mode,job=job,poses=len(rmsd),top1_rmsd_A=float(rmsd.iloc[0].heavy_atom_rmsd_A),
                best_rmsd_A=float(rmsd.heavy_atom_rmsd_A.min()),vina_score_kcal_mol=float(top.vina_score_kcal_mol),
                CNNscore=None if pd.isna(top.CNNscore) else float(top.CNNscore)))
    pd.DataFrame(gnina).to_csv(OUT/'gnina_comparison.csv',index=False)
    print(json.dumps(dict(haddock=rows,gnina=gnina),indent=2))



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
