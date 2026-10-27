"""Repack and minimize 1YCR sidechains under a fixed backbone for a small variant panel."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    import pyrosetta
    from pyrosetta import rosetta
    pyrosetta.init('-mute all -constant_seed -jran 20261020 -ex1 -ex2aro')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'data/peptide/1YCR_original.pdb'
    cleaned = output / '1YCR_protein.pdb'
    cleaned.write_text('\n'.join(line for line in source.read_text().splitlines()
                                  if line.startswith(('ATOM  ', 'TER', 'END'))) + '\n')
    original = pyrosetta.pose_from_pdb(str(cleaned))
    assert original.total_residue() == 98 and original.num_chains() == 2
    assert original.chain_sequence(2) == 'ETFSDLWKLLPEN'
    scorefxn = pyrosetta.create_score_function('ref2015')
    move_map = rosetta.core.kinematics.MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(True)
    move_map.set_jump(False)
    rows = []
    variants = [('WT', None, None), ('D5E', 21, 'GLU'), ('S4E', 20, 'GLU'), ('L9K', 25, 'LYS')]
    for name, pdb_residue, residue_name in variants:
        rosetta.numeric.random.rg().set_seed(20261020)
        pose = original.clone()
        position = pose.pdb_info().pdb2pose('B', pdb_residue) if pdb_residue else 0
        if pdb_residue:
            assert position > 0
            rosetta.protocols.simple_moves.MutateResidue(position, residue_name).apply(pose)
        before = scorefxn(pose)
        task = pyrosetta.standard_packer_task(pose)
        task.restrict_to_repacking()
        task.or_include_current(True)
        rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task).apply(pose)
        minimizer = rosetta.protocols.minimization_packing.MinMover()
        minimizer.movemap(move_map)
        minimizer.score_function(scorefxn)
        minimizer.min_type('lbfgs_armijo_nonmonotone')
        minimizer.tolerance(0.001)
        minimizer.apply(pose)
        after = scorefxn(pose)
        assert pose.num_jump() == 1
        interface = rosetta.protocols.analysis.InterfaceAnalyzerMover(1)
        interface.set_scorefunction(scorefxn)
        interface.set_pack_input(False)
        interface.set_pack_separated(True)
        interface.apply(pose)
        rows.append(dict(candidate=name, mutation_pdb=f'B{pdb_residue}{residue_name}' if pdb_residue else 'WT',
                         pose_position=position, peptide_sequence=pose.chain_sequence(2),
                         total_score_before_REU=before, total_score_after_REU=after,
                         interface_score_difference_REU=interface.get_interface_dG(),
                         interface_delta_sasa_A2=interface.get_interface_delta_sasa()))
        pose.dump_pdb(str(output / f'{name}.pdb'))
    with (output / 'scores.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    record = dict(completed=True, date='2026-10-20', pyrosetta_version=pyrosetta.version(),
                  scorefunction='ref2015', backbone_fixed=True, jumps_fixed=True,
                  sidechains_repacked='all protein residues', interface_jump=1,
                  pack_separated=True, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  candidates=rows, experimental_validation=False)
    (output / 'run_summary.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
