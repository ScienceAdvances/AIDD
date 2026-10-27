"""Audit actual filter flags rather than treating a ranked output as a pass."""
import json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'examples/debian13_20261023'
SOURCE=ROOT/'examples/debian13_20261023'

def main():
    d=pd.read_csv(SOURCE/'boltzgen_gpu/final_ranked_designs/final_designs_metrics_2.csv')
    cols=['id','final_rank','design_to_target_iptm','filter_rmsd','filter_rmsd_design','bb_target_aligned_rmsd_design','pass_filters','pass_filter_rmsd_filter','pass_filter_rmsd_design_filter']
    d[cols].to_csv(OUT/'boltzgen_candidate_summary.csv',index=False)
    r=dict(candidates=len(d), ranked=len(d),passed_all_filters=int(d.pass_filters.sum()),
           complex_filter_rmsd_range_A=[float(d.filter_rmsd.min()),float(d.filter_rmsd.max())],
           binder_only_rmsd_range_A=[float(d.filter_rmsd_design.min()),float(d.filter_rmsd_design.max())],
           esmfold2_executed=False, version='0.3.2',
           interpretation='Ranking penalizes failed filters but still fills the requested budget. Neither candidate passes all filters; binder folding alone does not establish a preserved complex pose.')
    (OUT/'boltzgen_analysis_summary.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))

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
