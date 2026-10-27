"""Use official RFD2 pipeline stages with explicitly selected native interpreters.

Only container launch commands are replaced. Model settings and scientific
implementations remain those in the explicitly configured RFdiffusion2 checkout.
"""
import builtins
import os
from pathlib import Path
import shlex
import sys
import runpy

ROOT = Path(__file__).resolve().parents[1]
RFD = Path(os.environ['RFD_REPO']).expanduser().resolve()
sys.path.insert(0, str(RFD))
from rf_diffusion.benchmark import mpnn_designs_v2, score_designs, add_metrics
from rf_diffusion.benchmark.util import slurm_tools

ligand_py = Path(os.environ['AIDD_LIGAND_PYTHON'])
chai_py = Path(os.environ['AIDD_CHAI_PYTHON'])
runtime = Path(os.environ['AIDD_RUNTIME'])
compat = ROOT / 'scripts/run_ligand_compat.py'
q = shlex.quote

mpnn_designs_v2.get_binary = lambda in_proc: f'{q(str(ligand_py))} {q(str(compat))}'


def translate(value):
    if isinstance(value, str):
        old = f'/usr/bin/apptainer run --nv {RFD}/rf_diffusion/exec/chai.sif {RFD}/rf_diffusion/benchmark/../../lib/chai/predict.py'
        value = value.replace(old, f'{q(str(chai_py))} {q(os.environ['AIDD_CHAI_PREDICT'])}')
        old = f'apptainer exec --env WANDB_MODE=offline {RFD}/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif python'
        value = value.replace(old, f'WANDB_MODE=offline {q(str(ligand_py))} {q(str(compat))}')
    return value


def native_print(*values, **kwargs):
    builtins.print(*(translate(value) for value in values), **kwargs)


score_designs.print = native_print
add_metrics.print = native_print
original_prune = slurm_tools.prune_jobs
slurm_tools.prune_jobs = lambda filename, checks: original_prune(
    filename, {translate(command): check for command, check in checks.items()})

if __name__ == '__main__':
    runpy.run_path(str(RFD / 'rf_diffusion/benchmark/pipeline.py'), run_name='__main__')
