"""Run the Boltz CLI with the 2.2.1 affinity-parser namespace bridge."""
from rdkit.Chem import AllChem, Descriptors

# Preserve the official molecular-weight implementation under its expected name.
AllChem.Descriptors = Descriptors

from boltz.main import cli

if __name__ == '__main__':
    cli()
