"""Run the Boltz 2.2.1 official YAML parser with trusted local CCD components."""
import argparse
import dataclasses
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ccd-zip', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from boltz.data import const
    from boltz.data.mol import load_canonicals
    from boltz.data.parse.yaml import parse_yaml
    from rdkit.Chem import AllChem, Descriptors
    # Boltz 2.2.1 accesses this legacy namespace during affinity parsing.
    AllChem.Descriptors = Descriptors
    args.output.mkdir(parents=True, exist_ok=False)
    moldir = args.output / 'canonical_molecules'
    moldir.mkdir()
    with zipfile.ZipFile(args.ccd_zip) as archive:
        entries = {Path(name).name: name for name in archive.namelist()}
        for token in const.canonical_tokens:
            name = f'{token}.pkl'
            (moldir / name).write_bytes(archive.read(entries[name]))
    ccd = load_canonicals(str(moldir))
    targets = []
    for filename in ['human_mdm2_p53_no_msa.yaml', 'official_affinity.yaml']:
        parsed = parse_yaml(ROOT / 'workflows/modern_tools/boltz' / filename, ccd, moldir, boltz2=True)
        targets.append(dict(input=filename, record=dataclasses.asdict(parsed.record),
                            atoms=len(parsed.structure.atoms), residues=len(parsed.structure.residues)))
    import rdkit
    result = dict(completed=True, boltz_package='2.2.1', rdkit_version=rdkit.__version__,
                  compatibility_bridge='AllChem.Descriptors = rdkit.Chem.Descriptors',
                  stage='official input parser', model_inference=False, targets=targets)
    (args.output / 'parser_summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
