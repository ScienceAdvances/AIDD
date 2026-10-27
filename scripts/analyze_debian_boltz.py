"""Compare new Boltz MDM2/p53 predictions to 1YCR and report affinity units."""
import csv
import json
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'examples/debian13_20261023'


SOURCE = ROOT / 'examples/debian13_20261023'

def chain(model, name):
    residues = [r for r in model[name] if r.id[0] == ' ' and 'CA' in r]
    return ''.join(seq1(r.resname) for r in residues), [r['CA'] for r in residues]


def main():
    ref = PDBParser(QUIET=True).get_structure('ref', ROOT / 'data/peptide/1YCR_original.pdb')[0]
    sa, ra = chain(ref, 'A'); sb, rb = chain(ref, 'B')
    rows = []
    for p in sorted((SOURCE / 'boltz_structure').rglob('*_model_*.cif')):
        model = MMCIFParser(QUIET=True).get_structure(p.stem, p)[0]
        a, ma = chain(model, 'A'); b, mb = chain(model, 'B')
        assert (a, b) == (sa, sb)
        align = Superimposer(); align.set_atoms(ra, ma); align.apply(model.get_atoms())
        rmsd = float(np.sqrt(np.mean(np.sum((np.array([x.coord for x in mb]) - np.array([x.coord for x in rb]))**2, axis=1))))
        confidence = json.loads((p.parent / ('confidence_' + p.stem + '.json')).read_text())
        rows.append(dict(model=p.stem, confidence_score=confidence.get('confidence_score'),
            iptm=confidence.get('iptm'), ptm=confidence.get('ptm'),
            mdm2_ca_rmsd_A=float(align.rms), p53_ca_rmsd_after_mdm2_alignment_A=rmsd))
    if rows:
        assert len(rows) == 5
        with (OUT / 'boltz_structure_comparison.csv').open('w', newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    affinity=[]
    for directory in ['boltz_affinity', 'boltz_affinity_no_msa']:
        for p in (SOURCE / directory).rglob('affinity_*.json'):
            record=json.loads(p.read_text())
            if 'affinity_pred_value' not in record: continue
            y=float(record['affinity_pred_value'])
            affinity.append(dict(source=str(p), input_msa='server' if directory=='boltz_affinity' else 'explicitly empty',
                predicted_log10_IC50_uM=y, predicted_IC50_uM=10**y, predicted_IC50_nM=1000*10**y,
                binding_probability=record.get('affinity_probability_binary'), experimental_validation=False))
    result=dict(structure_samples=len(rows), structure_rows=rows, affinity=affinity,
        experimental_validation=False, structure_reference_used_only_after_prediction=True)
    (OUT / 'boltz_analysis_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))



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
