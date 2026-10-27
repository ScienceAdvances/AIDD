"""Focused helper checks, using synthetic fixtures rather than docking results."""
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import csv
import io
import json
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem
import gnina_workflow as workflow


class WorkflowChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(self.molecule, randomSeed=17), 0)

    def tearDown(self):
        self.temporary.cleanup()

    def write_sdf(self, name, molecules):
        target = self.root / name
        writer = Chem.SDWriter(str(target))
        for molecule in molecules:
            writer.write(molecule)
        writer.close()
        return target

    def test_rmsd_keeps_receptor_frame(self):
        reference = self.write_sdf("reference.sdf", [self.molecule])
        moved = Chem.Mol(self.molecule)
        conformer = moved.GetConformer()
        for index in range(moved.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            conformer.SetAtomPosition(index, (point.x + 3.0, point.y, point.z))
        poses = self.write_sdf("poses.sdf", [self.molecule, moved])
        output = self.root / "rmsd.csv"
        with redirect_stdout(io.StringIO()):
            workflow.rmsd(SimpleNamespace(reference=reference, poses=poses, output=output))
        with output.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertAlmostEqual(float(rows[0]["heavy_atom_rmsd_A"]), 0.0, places=5)
        self.assertAlmostEqual(float(rows[1]["heavy_atom_rmsd_A"]), 3.0, places=4)
        self.assertEqual(rows[1]["within_2A"], "False")

    def test_rmsd_rejects_multiple_references(self):
        reference = self.write_sdf("reference.sdf", [self.molecule, self.molecule])
        poses = self.write_sdf("poses.sdf", [self.molecule])
        with self.assertRaisesRegex(ValueError, "exactly one"):
            workflow.rmsd(SimpleNamespace(reference=reference, poses=poses,
                                          output=self.root / "rmsd.csv"))

    def test_summary_preserves_all_poses_and_missing_scores(self):
        scored = Chem.Mol(self.molecule)
        scored.SetProp("minimizedAffinity", "-7.5")
        scored.SetProp("CNNscore", "0.8")
        self.write_sdf("fixture.sdf", [scored, self.molecule])
        manifest = [{"dry_run": False, "success": True,
                     "job": {"job_id": "fixture", "seed": 17,
                             "receptor": "synthetic_receptor.pdb", "ligand": "synthetic.sdf"}}]
        (self.root / "run_manifest.json").write_text(json.dumps(manifest))
        output = self.root / "scores.csv"
        with redirect_stdout(io.StringIO()):
            workflow.summarize(SimpleNamespace(run=self.root, output=output))
        with output.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["vina_score_kcal_mol"], "-7.5")
        self.assertEqual(rows[0]["CNNscore"], "0.8")
        self.assertEqual(rows[1]["vina_score_kcal_mol"], "")
        self.assertEqual(rows[1]["CNNscore"], "")
        self.assertEqual(rows[0]["CNNaffinity_pK"], "")
        scored.SetProp("CNNscore", "nan")
        with self.assertRaisesRegex(ValueError, "Non-finite"):
            workflow.property_value(scored, ["CNNscore"])

    def test_missing_output_is_failure_despite_zero_exit(self):
        reference = self.write_sdf("reference.sdf", [self.molecule])
        ligand = self.write_sdf("ligand.sdf", [self.molecule])
        receptor = self.root / "receptor.pdb"
        receptor.write_text("REMARK Synthetic execution fixture\n")
        binary = self.root / "fake_gnina"
        binary.write_text("#!/usr/bin/env python3\nimport sys\n"
                          "if '--version' in sys.argv:\n    print('synthetic-test-only')\n")
        binary.chmod(0o755)
        jobs = self.root / "jobs.tsv"
        with jobs.open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["job_id", "receptor", "ligand", "reference", "seed"])
            writer.writerow(["fixture", receptor, ligand, reference, 17])
        output = self.root / "execution"
        arguments = SimpleNamespace(cpu=1, exhaustiveness=1, num_modes=1,
                                    jobs=jobs, output=output, gnina=str(binary),
                                    dry_run=False, cnn_scoring="rescore", device=None)
        with self.assertRaisesRegex(RuntimeError, "1 jobs failed"):
            workflow.run(arguments)
        record = json.loads((output / "run_manifest.json").read_text())[0]
        self.assertEqual(record["returncode"], 0)
        self.assertFalse(record["success"])
        self.assertIn("error", record)


if __name__ == "__main__":
    unittest.main()
