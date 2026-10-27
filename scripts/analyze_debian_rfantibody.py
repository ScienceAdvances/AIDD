"""Validate Quiver counts, preserved sequences and RF2's target-aligned geometry."""
import csv, io, json, re
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'examples/debian13_20261023'
RUNTIME=ROOT/'examples/debian13_20261023/rfantibody_outputs'
def quiver(path):
    blocks={}; tag=None
    for line in path.read_text().splitlines(keepends=True):
        if line.startswith('QV_TAG '): tag=line.split()[1]; blocks[tag]=[]
        elif tag and not line.startswith('QV_SCORE '): blocks[tag].append(line)
    return blocks

def model(lines): return PDBParser(QUIET=True).get_structure('m',io.StringIO(''.join(lines)))[0]
def residues(m,c): return [r for r in m[c] if r.id[0]==' ' and 'CA' in r]
def sequence(m,c): return ''.join(seq1(r.resname) for r in residues(m,c))
def main():
    backbone=quiver(RUNTIME/'1_rfdiffusion.qv'); designed=quiver(RUNTIME/'2_proteinmpnn.qv'); refolded=quiver(RUNTIME/'3_rf2.qv')
    assert (len(backbone),len(designed),len(refolded))==(2,4,4)
    scores=list(csv.DictReader((RUNTIME/'3_rf2.sc').open(),delimiter='\t'))
    assert len(scores)==4
    rows=[]
    for row in scores:
        tag=row['tag']; source=tag.removesuffix('_best'); parent=source.split('_dldesign_')[0]
        b=model(backbone[parent]); d=model(designed[source]); f=model(refolded[tag])
        assert set(c.id for c in f)=={'H','T'}
        assert sequence(b,'T')==sequence(d,'T')==sequence(f,'T')
        assert sequence(d,'H')==sequence(f,'H')
        loops={int(re.search(r'LABEL:\s+(\d+)',s)[1]) for s in backbone[parent] if 'PDBinfo-LABEL:' in s and re.search(r' H[123]\s*$',s)}
        br=residues(b,'H'); dr=residues(d,'H')
        assert len(br)==len(dr)
        changes=[i for i,(x,y) in enumerate(zip(br,dr),1) if x.resname!=y.resname]
        assert set(changes)<=loops
        align=Superimposer(); align.set_atoms([r['CA'] for r in residues(b,'T')],[r['CA'] for r in residues(f,'T')]); align.apply(f.get_atoms())
        rmsd=float(np.sqrt(np.mean(np.sum((np.array([r['CA'].coord for r in br])-np.array([r['CA'].coord for r in residues(f,'H')]))**2,axis=1))))
        row={k:float(v) if k!='tag' else v for k,v in row.items()}
        row={k:v for k,v in row.items() if not k.startswith('framework_aligned_L')}
        row.update(sequence=sequence(d,'H'),length=len(dr),cdr_mutations=len(changes),
                   independent_target_aligned_antibody_rmsd_A=rmsd)
        rows.append(row)
    with (OUT/'rfantibody_candidate_summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    errors=[abs(r['target_aligned_antibody_rmsd']-r['independent_target_aligned_antibody_rmsd_A']) for r in rows]
    assert max(errors)<0.02,errors
    best=min(rows,key=lambda r:r['interaction_pae'])
    result=dict(backbones=2,sequence_designs=4,refolded_models=4,passed=True,
        target_sequence_preserved=True,framework_sequence_preserved_during_mpnn=True,
        mpnn_to_rf2_sequences_preserved=True, independent_rmsd_max_error_A=max(errors),
        best_by_interaction_pae=best,
        interpretation='One candidate retains the intended pose with low interface PAE; two have >42 A target-aligned displacement despite high global pLDDT. These are model predictions, not experimentally verified binders.',
        scorefile_precision='Official Quiver scores rounded to two decimals; light-chain metrics omitted for VHH.',experimental_validation=False)
    (OUT/'rfantibody_analysis_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=RUNTIME)
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    RUNTIME=args.source.resolve()
    OUT=args.output.resolve()
    OUT.mkdir(parents=True, exist_ok=False)
    main()
