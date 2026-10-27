"""Compare AF3 samples with 1YCR after aligning the MDM2 C-alpha atoms."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]


def ca_atoms(model, chain):
    residues = [residue for residue in model[chain] if residue.id[0] == ' ' and 'CA' in residue]
    return ''.join(seq1(r.resname) for r in residues), [r['CA'] for r in residues]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=args.overwrite)
    ranking_files = list(args.run.glob('*ranking_scores.csv'))
    exact_ranking = {}
    if len(ranking_files) == 1:
        with ranking_files[0].open() as stream:
            exact_ranking = {f"seed-{row['seed']}_sample-{row['sample']}": float(row['ranking_score'])
                             for row in csv.DictReader(stream)}
    reference = PDBParser(QUIET=True).get_structure('reference', ROOT / 'data/peptide/1YCR_original.pdb')[0]
    ref_seq_a, ref_a = ca_atoms(reference, 'A')
    ref_seq_b, ref_b = ca_atoms(reference, 'B')
    rows = []
    for folder in sorted(args.run.glob('seed-*_sample-*')):
        model_path = next(folder.glob('*_model.cif'))
        model = MMCIFParser(QUIET=True).get_structure('prediction', model_path)[0]
        sequence_a, moving_a = ca_atoms(model, 'A')
        sequence_b, moving_b = ca_atoms(model, 'B')
        assert (sequence_a, sequence_b) == (ref_seq_a, ref_seq_b)
        align = Superimposer()
        align.set_atoms(ref_a, moving_a)
        align.apply(model.get_atoms())
        ligand_rmsd = np.sqrt(np.mean(np.sum((np.array([a.coord for a in moving_b]) - np.array([a.coord for a in ref_b]))**2, axis=1)))
        confidence = json.loads(next(folder.glob('*_summary_confidences.json')).read_text())
        details_path = next(p for p in folder.glob('*_confidences.json') if not p.name.endswith('_summary_confidences.json'))
        details = json.loads(details_path.read_text())
        atom_plddt = np.array(details['atom_plddts'])
        atom_chains = np.array(details['atom_chain_ids'])
        rows.append(dict(sample=folder.name, ptm=confidence['ptm'], iptm=confidence['iptm'],
                         ranking_score=exact_ranking.get(folder.name, confidence['ranking_score']), has_clash=confidence['has_clash'],
                         mdm2_ca_rmsd_A=float(align.rms), p53_ca_rmsd_after_mdm2_alignment_A=float(ligand_rmsd),
                         mdm2_atom_mean_plddt=float(atom_plddt[atom_chains == 'A'].mean()),
                         p53_atom_mean_plddt=float(atom_plddt[atom_chains == 'B'].mean())))
    assert rows
    with (args.output / 'sample_comparison.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    best = max(rows, key=lambda row: row['ranking_score'])
    summary = dict(completed=True, reference='human 1YCR', reference_coordinates_used_as_prediction_input=False,
                   msa='explicitly empty', samples=len(rows), best_by_af3_ranking=best,
                   metric='p53 CA RMSD after least-squares alignment on 85 MDM2 CA atoms', rows=rows)
    (args.output / 'analysis_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    best_dir = args.run / best['sample']
    details_path = next(p for p in best_dir.glob('*_confidences.json') if not p.name.endswith('_summary_confidences.json'))
    details = json.loads(details_path.read_text())
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=True)
    im = axes[0].imshow(details['pae'], vmin=0, vmax=31.75, cmap='viridis_r')
    axes[0].axhline(84.5, color='white', linewidth=0.8); axes[0].axvline(84.5, color='white', linewidth=0.8)
    axes[0].set(title='Best-ranked sample: predicted aligned error', xlabel='Token', ylabel='Token')
    fig.colorbar(im, ax=axes[0], label='PAE (A)')
    axes[1].bar(range(len(rows)), [row['p53_ca_rmsd_after_mdm2_alignment_A'] for row in rows], color='#327a88')
    axes[1].set(title='p53 displacement after MDM2 alignment', xlabel='Sample index', ylabel='CA RMSD (A)')
    fig.savefig(args.output / 'af3_pae_and_reference_comparison.png', dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
