"""Check official motif atoms, antibody labels, prediction inputs, and VHH records."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'workflows/modern_tools'


def atoms(path):
    return [dict(record=line[:6].strip(), atom=line[12:16].strip(), residue=line[17:20].strip(),
                 chain=line[21], number=int(line[22:26]))
            for line in path.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    rfd_atoms = atoms(INPUT / 'rfdiffusion2/M0584_1ldm.pdb')
    motif = {106: ['NE', 'CD', 'CZ'], 166: ['OD1', 'CG'], 169: ['NH2', 'CZ'], 193: ['NE2', 'CD2', 'CE1']}
    for residue, names in motif.items():
        available = {a['atom'] for a in rfd_atoms if a['chain'] == 'A' and a['number'] == residue}
        assert set(names) <= available, (residue, names, available)
    ligand_atoms = {name: {a['atom'] for a in rfd_atoms if a['residue'] == name} for name in ['NAD', 'OXM']}
    assert {'O7N', 'C7N', 'C3N', 'N7N', 'C2N', 'C4N', 'N1N', 'C5N', 'C1D'} <= ligand_atoms['NAD']
    assert {'O3', 'C2', 'C1', 'O2', 'N1'} <= ligand_atoms['OXM']
    target = atoms(INPUT / 'rfantibody/flu_HA.pdb')
    framework = atoms(INPUT / 'rfantibody/h-NbBCII10.pdb')
    for residue in [146, 170, 177]:
        assert any(a['chain'] == 'B' and a['number'] == residue and a['atom'] == 'CA' for a in target)
    labels = Counter()
    for line in (INPUT / 'rfantibody/h-NbBCII10.pdb').read_text().splitlines():
        if line.startswith('REMARK PDBinfo-LABEL:'):
            labels[line.split()[-1]] += 1
    assert set(labels) == {'H1', 'H2', 'H3'}
    af3 = json.loads((INPUT / 'alphafold3/human_mdm2_p53_no_msa.json').read_text())
    lengths = {s['protein']['id']: len(s['protein']['sequence']) for s in af3['sequences']}
    assert lengths == {'A': 85, 'B': 13}
    results = dict(date=__import__('datetime').date.today().isoformat(), checks_passed=True,
                   rfdiffusion2=dict(motif_residues=4, motif_atoms=sum(map(len, motif.values())),
                                    ligand_atom_counts={k: len(v) for k, v in ligand_atoms.items()}, design_executed=False),
                   rfantibody=dict(target_chain_residues=dict(Counter(a['chain'] for a in target if a['atom'] == 'CA')),
                                  framework_chain_residues=dict(Counter(a['chain'] for a in framework if a['atom'] == 'CA')),
                                  cdr_label_counts=dict(labels), hotspots_checked=['B146', 'B170', 'B177'], design_executed=False),
                   alphafold3_input=dict(chain_lengths=lengths, msa='explicitly empty', templates=0))
    (args.output / 'input_checks.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
