"""Build a traceable study brief from saved molecular-design case outputs."""
import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequence-run', type=Path, default=ROOT / 'examples/debian13_20261023/sequence_debian_only')
    parser.add_argument('--cart-run', type=Path, default=ROOT / 'examples/debian13_20261023/cart')
    parser.add_argument('--md-analysis', type=Path, default=ROOT / 'examples/debian13_20261023/md_analysis/analysis_summary.json')
    parser.add_argument('--gnina-run', type=Path, default=ROOT / 'examples/debian13_20261023/gnina_rescore_gpu')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--haddock-table', type=Path, default=ROOT / 'examples/debian13_20261023/haddock_analysis/full_final_models.csv')
    parser.add_argument('--haddock-label', default='HADDOCK Debian full')
    args = parser.parse_args()
    sequence_run, cart_run = args.sequence_run.resolve(), args.cart_run.resolve()
    summary = json.loads((sequence_run / 'run_summary.json').read_text())
    if not summary['completed']:
        raise ValueError('The sequence run is incomplete')
    candidates = read_csv(sequence_run / 'candidates.csv')
    shortlist = read_csv(sequence_run / 'shortlist.csv')
    cart = read_csv(cart_run / 'panel_summary.csv')
    md = json.loads(args.md_analysis.read_text())
    gnina = json.loads((args.gnina_run / 'run_manifest.json').read_text())
    haddock_path = args.haddock_table.resolve()
    haddock = sorted(read_csv(haddock_path), key=lambda row: float(row['score']))
    best = max(haddock, key=lambda row: float(row['dockq']))
    evidence = []

    def add(case, metric, value, unit, source, interpretation, next_step):
        evidence.append(dict(case=case, metric=metric, value=value, unit=unit,
                             source=str(source), interpretation=interpretation, next_step=next_step))

    add(args.haddock_label, 'top_score_dockq', haddock[0]['dockq'], 'dimensionless',
        haddock_path, 'Compare the best scoring pose with the homologous reference using DockQ.',
        'Inspect scoring, interfaces and reference mapping.')
    add(args.haddock_label, 'best_output_dockq', best['dockq'], 'dimensionless',
        haddock_path, f'Best reference pose ranks {haddock.index(best) + 1} by score.',
        'Separate sampling coverage from selection performance.')
    wt = next(row for row in candidates if row['candidate_id'] == 'WT')
    for row in shortlist:
        add('Sequence optimization', row['candidate_id'] + '_nll', row['mpnn_nll_mean'],
            'mean negative natural log probability', sequence_run / 'candidates.csv',
            'Fixed-backbone sequence compatibility under a common scoring protocol.',
            'Recheck structure, preparation properties and binding.')
    for row in cart:
        if row['panel'] in ('4H', 'S7B', 'S7I_5:1') and row['group'] in ('[B]-L3', '[B]-H2L3wt'):
            unit = {'4H': 'relative IFNg', 'S7B': 'MFI', 'S7I_5:1': 'relative viability percent'}[row['panel']]
            add('CAR-T', row['panel'] + '_' + row['group'], row['mean'], unit,
                cart_run / 'panel_summary.csv', 'Panel-specific descriptive mean.',
                'Assess expression, binding, background and function together.')
    for row in md:
        add('Factor Xa MD', row['source'] + '_rmsd_to_em', row['rmsd_to_em_mean_nm'], 'nm',
            args.md_analysis.resolve(), f"Observed over {row['end_ps'] - row['start_ps']:.1f} ps.",
            'Review equilibration, sampling and independent repeats.')
    dry = all(row['dry_run'] for row in gnina)
    add('GNINA', 'workflow_state', 'dry_run' if dry else 'engine_run', 'state',
        args.gnina_run.resolve() / 'run_manifest.json',
        'Commands and inputs checked.' if dry else 'Inspect each recorded success state and pose output.',
        'Run the engine and evaluate pose RMSD.' if dry else 'Evaluate top-1, best-output RMSD and seed variation.')
    hypotheses = {
        'WT': 'Baseline for preparation and binding comparisons.',
        'single_D5E': 'Conservative acidic substitution at one position.',
        'single_S4E': 'A local acidic substitution with a lower sequence hydrophobicity proxy.',
        'single_L9K': 'Lower hydrophobicity proxy with a possible interface-contact tradeoff.',
        'single_L9R': 'Alternative charged substitution at the same local position.',
    }
    plan = [dict(candidate_id=row['candidate_id'], sequence=row['sequence'],
                 delta_nll_vs_wt=float(row['mpnn_nll_mean']) - float(wt['mpnn_nll_mean']),
                 hypothesis=hypotheses.get(row['candidate_id'], 'Inspect the actual substitution and local environment.'),
                 structure_check='Re-predict complex; inspect hotspot geometry and clashes.',
                 preparation_check='Confirm identity, purity, solubility and chemical state.',
                 binding_check='Compare with WT under matched binding-assay conditions.') for row in shortlist]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / 'evidence.csv', evidence)
    write_csv(output / 'candidate_validation_plan.csv', plan)
    lines = ['# Case study brief', '',
             '## Research question', '',
             'How can structural and sequence models prioritize variants for further validation?', '',
             '## Evidence', '',
             '| Case | Metric | Value | Unit |', '|---|---|---|---|']
    lines.extend(f"| {row['case']} | {row['metric']} | {row['value']} | {row['unit']} |" for row in evidence)
    lines += ['', '## Decision and next work', '',
              f"The sequence run contains {len(candidates)} candidates and {len(shortlist)} nominated sequences, including WT.",
              'Use candidate_validation_plan.csv to connect each nomination to a testable hypothesis.',
              'Use evidence.csv to recover the measurement definition and source of each claim.',
              'Describe Factor Xa MD and CAR-T reanalysis as separate case studies with their own biological contexts.',
              'GNINA execution state follows its run manifest.', '']
    (output / 'brief.md').write_text('\n'.join(lines))
    (output / 'build_summary.json').write_text(json.dumps(dict(
        completed=True, evidence_rows=len(evidence), candidate_plan_rows=len(plan),
        sequence_run=str(sequence_run), cart_run=str(cart_run),
        md_analysis=str(args.md_analysis.resolve()), gnina_run=str(args.gnina_run.resolve())), indent=2) + '\n')
    print(json.dumps(dict(output=str(output), evidence_rows=len(evidence), candidate_plan_rows=len(plan))))


if __name__ == '__main__':
    main()
