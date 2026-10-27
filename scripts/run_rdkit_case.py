"""Prepare an imatinib conformer ensemble and characterize the p53 variant panel."""
import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'workflows/gnina/example_1iep/reference.sdf')
    parser.add_argument('--num-confs', type=int, default=20)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.num_confs < 1: parser.error("--num-confs must be positive")
    import rdkit
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator, rdMolDescriptors
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    supplier = Chem.SDMolSupplier(str(args.input), removeHs=False)
    reference = next(iter(supplier))
    if reference is None:
        raise ValueError('The reference SDF is invalid')
    molecule = Chem.RemoveHs(reference)
    canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
    generated = Chem.AddHs(molecule)
    generated.RemoveAllConformers()
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20261020
    parameters.numThreads = 1
    parameters.pruneRmsThresh = 0.5
    conformers = list(AllChem.EmbedMultipleConfs(generated, numConfs=args.num_confs, params=parameters))
    if not conformers or not AllChem.MMFFHasAllMoleculeParams(generated):
        raise RuntimeError('Conformer generation or MMFF parameterization failed')
    energies = AllChem.MMFFOptimizeMoleculeConfs(generated, numThreads=1, maxIters=2000)
    rows = []
    with Chem.SDWriter(str(output / 'imatinib_conformers.sdf')) as writer:
        for conformer, (status, energy) in zip(conformers, energies):
            generated.SetProp('_Name', f'imatinib_conformer_{conformer}')
            generated.SetIntProp('mmff_status', status)
            generated.SetDoubleProp('mmff_energy_kcal_mol', energy)
            writer.write(generated, confId=conformer)
            rows.append(dict(conformer_id=conformer, mmff_status=status, mmff_energy_kcal_mol=energy))
    with (output / 'conformer_energies.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprint = generator.GetFingerprint(molecule)
    roundtrip = Chem.MolFromSmiles(canonical)
    similarity = DataStructs.TanimotoSimilarity(fingerprint, generator.GetFingerprint(roundtrip))
    panel = []
    for name, sequence in [('WT', 'ETFSDLWKLLPEN'), ('D5E', 'ETFSELWKLLPEN'),
                           ('S4E', 'ETFEDLWKLLPEN'), ('L9K', 'ETFSDLWKKLPEN')]:
        peptide = Chem.MolFromFASTA(sequence)
        if peptide is None:
            raise ValueError(f'Invalid peptide: {name}')
        panel.append(dict(candidate=name, sequence=sequence, formula=rdMolDescriptors.CalcMolFormula(peptide),
                          molecular_weight_Da=Descriptors.MolWt(peptide),
                          formal_charge_as_constructed=Chem.GetFormalCharge(peptide),
                          heavy_atoms=peptide.GetNumHeavyAtoms(),
                          isomeric_smiles=Chem.MolToSmiles(peptide, isomericSmiles=True)))
    with (output / 'peptide_properties.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(panel[0]));writer.writeheader();writer.writerows(panel)
    summary = dict(completed=True, date=__import__('datetime').datetime.now().astimezone().isoformat(), rdkit_version=rdkit.__version__,
                   input=str(args.input.resolve()), canonical_smiles=canonical,
                   heavy_atoms=molecule.GetNumHeavyAtoms(), formal_charge=Chem.GetFormalCharge(molecule),
                   formula=rdMolDescriptors.CalcMolFormula(molecule), molecular_weight_Da=Descriptors.MolWt(molecule),
                   logP=Descriptors.MolLogP(molecule), tpsa_A2=rdMolDescriptors.CalcTPSA(molecule),
                   hbd=rdMolDescriptors.CalcNumHBD(molecule), hba=rdMolDescriptors.CalcNumHBA(molecule),
                   requested_conformers=args.num_confs, generated_conformers=len(conformers),
                   converged_conformers=sum(row['mmff_status'] == 0 for row in rows),
                   canonical_roundtrip_tanimoto=similarity, peptide_panel=panel,
                   ph_state_enumerated=False, protein_binding_evaluated=False)
    (output / 'run_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
