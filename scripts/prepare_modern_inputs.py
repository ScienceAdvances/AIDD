"""Assemble small, versioned teaching inputs from checked official sources."""
import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'workflows/modern_tools'
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))


def pdb_residues(path):
    residues = {}
    for line in path.read_text().splitlines():
        if line.startswith('ATOM  ') and line[12:16].strip() == 'CA':
            residues[(line[21], int(line[22:26]), line[26])] = AA[line[17:20]]
    return residues


def main():
    manifest = []
    previous_path = DEST / 'input_manifest.json'
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else {}
    expected = {entry['file']: entry['sha256'] for entry in previous.get('official_inputs', [])}
    copies = [
        ('rfdiffusion2/M0584_1ldm.pdb', 'rfdiffusion2/M0584_1ldm.pdb', 'https://raw.githubusercontent.com/RosettaCommons/RFdiffusion2/d365cbf4db3958814a9f8e4f6f94fa309dfebc2b/rf_diffusion/benchmark/input/mcsa_41/M0584_1ldm.pdb'),
        ('rfdiffusion2/checked_docs/rf_diffusion/benchmark/open_source_demo.json', 'rfdiffusion2/open_source_demo.json', 'https://raw.githubusercontent.com/RosettaCommons/RFdiffusion2/d365cbf4db3958814a9f8e4f6f94fa309dfebc2b/rf_diffusion/benchmark/open_source_demo.json'),
        ('rfantibody/checked_docs/scripts/examples/example_inputs/flu_HA.pdb', 'rfantibody/flu_HA.pdb', 'https://raw.githubusercontent.com/RosettaCommons/RFantibody/3e5a5bb1cb98befde9fc049ae3ff58cabb491f67/scripts/examples/example_inputs/flu_HA.pdb'),
        ('rfantibody/checked_docs/scripts/examples/example_inputs/h-NbBCII10.pdb', 'rfantibody/h-NbBCII10.pdb', 'https://raw.githubusercontent.com/RosettaCommons/RFantibody/3e5a5bb1cb98befde9fc049ae3ff58cabb491f67/scripts/examples/example_inputs/h-NbBCII10.pdb'),
        ('boltz/checked_docs/examples/affinity.yaml', 'boltz/official_affinity.yaml', 'https://raw.githubusercontent.com/jwohlwend/boltz/v2.2.1/examples/affinity.yaml'),
        ('', 'boltzgen/9bkq-assembly2.cif', 'https://raw.githubusercontent.com/HannesStark/boltzgen/a07fea1153111ef453aaea7a15697ee05de9291f/example/nanobody/9bkq-assembly2.cif'),
        ('', 'boltzgen/7xl0.cif', 'https://raw.githubusercontent.com/HannesStark/boltzgen/a07fea1153111ef453aaea7a15697ee05de9291f/example/nanobody_scaffolds/7xl0.cif'),
        ('', 'boltzgen/7xl0_official.yaml', 'https://raw.githubusercontent.com/HannesStark/boltzgen/a07fea1153111ef453aaea7a15697ee05de9291f/example/nanobody_scaffolds/7xl0.yaml'),
    ]
    for source, target, url in copies:
        path = DEST / target
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            with urllib.request.urlopen(url, timeout=60) as response:
                path.write_bytes(response.read())
        relative = str(path.relative_to(ROOT))
        actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if relative in expected and actual_hash != expected[relative]:
            raise ValueError(f'Official input changed: {relative}; save edits as a derived file')
        manifest.append(dict(file=relative, source_url=url, sha256=actual_hash))
    residues = pdb_residues(ROOT / 'data/peptide/1YCR_original.pdb')
    chains = {chain: ''.join(aa for (c, _, _), aa in residues.items() if c == chain) for chain in ['A', 'B']}
    assert len(chains['A']) == 85 and chains['B'] == 'ETFSDLWKLLPEN'
    af3 = dict(name='human_mdm2_p53_workflow', modelSeeds=[20261020],
               sequences=[dict(protein=dict(id=chain, sequence=sequence, unpairedMsa='', pairedMsa='', templates=[])) for chain, sequence in chains.items()],
               dialect='alphafold3', version=3)
    path = DEST / 'alphafold3/human_mdm2_p53_no_msa.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(af3, indent=2) + '\n')
    (DEST / 'boltz/human_mdm2_p53_no_msa.yaml').write_text('version: 1\nsequences:\n' + ''.join(f'  - protein:\n      id: {chain}\n      sequence: {sequence}\n      msa: empty\n' for chain, sequence in chains.items()))
    scaffold = (DEST / 'boltzgen/7xl0_official.yaml').read_text()
    scaffold = scaffold.replace('      id: A\ndesign:', '      id: A\n      res_index: 1..121\ndesign:', 1)
    assert 'res_index: 1..121' in scaffold
    (DEST / 'boltzgen/7xl0_resolved.yaml').write_text(scaffold)
    (DEST / 'boltzgen/penguinpox.yaml').write_text('entities:\n  - file:\n      path: 9bkq-assembly2.cif\n      include:\n        - chain:\n            id: B\n  - file:\n      path: 7xl0_resolved.yaml\n')
    derived = []
    for name in ['alphafold3/human_mdm2_p53_no_msa.json', 'boltz/human_mdm2_p53_no_msa.yaml',
                 'boltzgen/7xl0_resolved.yaml', 'boltzgen/penguinpox.yaml']:
        path = DEST / name
        derived.append(dict(file=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (DEST / 'input_manifest.json').write_text(json.dumps(dict(date='2026-10-20', official_inputs=manifest,
                                                           derived_inputs=derived, scaffold_edit='include chain A residues 1..121',
                                                           human_case=dict(pdb='1YCR', chains=chains, msa='explicitly empty for workflow exercise', experimental_coordinates_used_as_prediction_input=False)), indent=2) + '\n')
    print(json.dumps(dict(official_files=len(copies), human_chain_lengths={c: len(s) for c, s in chains.items()})))


if __name__ == '__main__':
    main()
