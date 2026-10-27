"""Preserve removed NumPy aliases for the bundled historical OpenFold code."""
import runpy
import sys
import numpy as np

np.int = int
np.bool = bool
script, *args = sys.argv[1:]
sys.argv = [script, *args]
sys.path.insert(0, str(__import__('pathlib').Path(script).parent))
runpy.run_path(script, run_name='__main__')
