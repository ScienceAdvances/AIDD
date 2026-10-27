"""Preview or execute a chapter using reader-supplied executable/repository paths."""
import argparse
import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]

def plan(chapter,c,out):
    tasks=[];py=c['analysis_python'];inp=ROOT/'workflows/modern_tools'
    def add(name,argv,cwd=ROOT,expected=(),env=None):tasks.append(dict(name=name,argv=list(map(str,argv)),cwd=str(cwd),expected=list(map(str,expected)),env=env or {}))
    def script(name,file,*args,python=py,expected=()):add(name,[python,ROOT/'scripts'/file,*args],expected=expected)
    def require(key):
        if not c.get(key):raise ValueError(f'Set {key} in config/tools.json; see docs/ENVIRONMENTS.md')
        return Path(c[key]).expanduser().resolve()
    if chapter=='01':
        work=out/'practice'
        script('prepare','prepare_haddock_practice.py','--output',work,expected=[work/'workflows/smoke.cfg'])
        add('haddock_smoke',[c['haddock3'],'workflows/smoke.cfg'],cwd=work,expected=[work/'runs/smoke/13_caprieval/capri_ss.tsv'])
    elif chapter=='02':
        script('design','run_proteinmpnn.py','--output',out/'design','--device','cuda',expected=[out/'design/candidate_summary.csv'])
        script('audit','analyze_proteinmpnn.py','--run',out/'design',expected=[out/'design/probability_audit.json'])
    elif chapter=='03':script('md','run_md.py','--gmx',c['gmx'],'--output',out/'md','--full','--gpu',expected=[out/'md/md.xtc',out/'md/run_summary.json'])
    elif chapter=='04':
        script('extract','extract_cart_data.py','--output',out/'cart',expected=[out/'cart/source_observations.csv'])
        script('plot','plot_cart.py','--input',out/'cart/source_observations.csv','--output',out/'cart',expected=[out/'cart/mutation_tradeoffs.png'])
    elif chapter=='06':
        script('score','optimize_sequence.py','--output',out/'scoring','--device','cuda','--repeats','8','--batch-size','16',expected=[out/'scoring/candidates.csv'])
        script('audit','verify_sequence_optimization.py','--run',out/'scoring',expected=[out/'scoring/verification.json'])
    elif chapter=='07':
        workflow=ROOT/'workflows/gnina/gnina_workflow.py'
        for mode in ['rescore','none']:
            dest=out/mode
            add(mode,[py,workflow,'run','--jobs',ROOT/'workflows/gnina/jobs_1iep.tsv','--output',dest,'--gnina',c['gnina'],'--device','0','--cpu','4','--exhaustiveness','32','--num-modes','20','--cnn-scoring',mode],expected=[dest/'run_manifest.json',*[dest/f'abl_imatinib_seed{s}.sdf' for s in [11,23,47]]])
            add(mode+'_summarize',[py,workflow,'summarize','--run',dest,'--output',dest/'scores.csv'],expected=[dest/'scores.csv'])
            for seed in [11,23,47]:add(f'{mode}_{seed}_rmsd',[py,workflow,'rmsd','--reference',ROOT/'workflows/gnina/example_1iep/reference.sdf','--poses',dest/f'abl_imatinib_seed{seed}.sdf','--output',dest/f'abl_imatinib_seed{seed}_rmsd.csv'])
    elif chapter=='09':script('rdkit','run_rdkit_case.py','--output',out/'rdkit',expected=[out/'rdkit/imatinib_conformers.sdf'])
    elif chapter=='10':script('pyrosetta','run_pyrosetta_case.py','--output',out/'pyrosetta',python=c['pyrosetta_python'],expected=[out/'pyrosetta/scores.csv',out/'pyrosetta/WT.pdb'])
    elif chapter=='11':
        rfd=require('rfd_repo');runtime=out/'runtime'
        env={'PYTHONPATH':str(rfd),'RFD_REPO':str(rfd),'AIDD_LIGAND_PYTHON':c['ligand_python'],'AIDD_CHAI_PYTHON':c['chai_python'],'AIDD_CHAI_PREDICT':str(require('chai_predict')),'AIDD_RUNTIME':str(runtime),'PATH':str(Path(c['rfd_python']).parent)+os.pathsep+os.environ.get('PATH','')}
        add('backbone',[c['rfd_python'],rfd/'rf_diffusion/benchmark/pipeline.py','--config-name=open_source_demo','sweep.benchmarks=active_site_unindexed_atomic_partial_ligand',f'outdir={runtime}/rfdiffusion2','sweep.retries=0'],cwd=runtime,env=env)
        add('ligandmpnn_chai',[c['ligand_python'],ROOT/'scripts/run_rfd_native_pipeline.py','--config-name=open_source_demo',f'outdir={runtime}/rfdiffusion2','start_step=mpnn','stop_step=end','score.chai1.omit_esm_embeddings=True'],cwd=runtime,env=env,expected=[runtime/'rfdiffusion2/compiled_metrics.csv'])
        script('analyze','analyze_debian_rfd.py','--source',runtime/'rfdiffusion2/compiled_metrics.csv','--output',out/'analysis',expected=[out/'analysis/rfdiffusion2_candidate_summary.csv'])
    elif chapter=='12':
        repo=require('rfantibody_repo');binary=require('rfantibody_bin');runtime=out/'runtime';dest=runtime/'rfantibody';python=c['rfantibody_python']
        env={'RFA_WEIGHTS_DIR':str(repo/'weights'),'CUBLAS_WORKSPACE_CONFIG':':4096:8','PATH':str(binary)+os.pathsep+os.environ.get('PATH','')}
        add('backbone',[binary/'rfdiffusion','--target',inp/'rfantibody/flu_HA.pdb','--framework',inp/'rfantibody/h-NbBCII10.pdb','--output-quiver',dest/'1_rfdiffusion.qv','--num-designs','2','--design-loops','H1:7,H2:6,H3:5-13','--hotspots','B146,B170,B177','--diffuser-t','50','--deterministic','--extra',f'hydra.run.dir={runtime}/hydra_rfantibody','--extra','hydra.job.chdir=True'],cwd=runtime,env=env,expected=[dest/'1_rfdiffusion.qv'])
        add('sequence',[python,repo/'scripts/proteinmpnn_interface_design.py','-quiver',dest/'1_rfdiffusion.qv','-outquiver',dest/'2_proteinmpnn.qv','-checkpoint_path',repo/'weights/ProteinMPNN_v48_noise_0.2.pt','-checkpoint_name',dest/'mpnn.check.point','-seqs_per_struct','2','-temperature','0.2','-debug','-deterministic'],cwd=runtime,env=env,expected=[dest/'2_proteinmpnn.qv'])
        add('rf2',[python,repo/'scripts/rf2_predict.py',f'input.quiver={dest}/2_proteinmpnn.qv',f'output.quiver={dest}/3_rf2.qv','inference.num_recycles=10','inference.cautious=True','inference.hotspot_show_proportion=0.0',f'model.model_weights={repo}/weights/RF2_ab.pt','+inference.seed=42',f'hydra.run.dir={runtime}/hydra_rf2','hydra.job.chdir=True'],cwd=runtime,env=env,expected=[dest/'3_rf2.qv'])
        add('extract',[binary/'qvextract',dest/'3_rf2.qv','-o',dest/'final_pdbs'],env=env,expected=[dest/'final_pdbs'])
        add('scores',[binary/'qvscorefile',dest/'3_rf2.qv'],env=env,expected=[dest/'3_rf2.sc'])
        script('analyze','analyze_debian_rfantibody.py','--source',dest,'--output',out/'analysis')
    elif chapter=='13':
        repo=require('af3_repo');params=require('af3_params');dest=out/'prediction'
        add('af3',[c['af3_python'],repo/'run_alphafold.py',f'--json_path={inp}/alphafold3/human_mdm2_p53_no_msa.json',f'--model_dir={params}',f'--output_dir={dest}','--run_data_pipeline=false','--flash_attention_implementation=xla','--buckets=128','--num_diffusion_samples=5'],expected=[dest/'human_mdm2_p53_workflow/human_mdm2_p53_workflow_ranking_scores.csv'])
        script('analyze','analyze_af3_case.py','--run',dest/'human_mdm2_p53_workflow','--output',out/'analysis')
    elif chapter=='14':
        common=['--model','boltz2','--cache',ROOT/'downloads/boltz','--accelerator','gpu','--devices','1','--seed','20261020','--num_workers','0','--max_parallel_samples','1','--recycling_steps','3','--sampling_steps','200','--diffusion_samples','5','--no_kernels','--write_full_pae','--write_full_pde']
        for label,name,extra in [('structure','human_mdm2_p53_no_msa',[]),('affinity','official_affinity',['--use_msa_server'])]:
            dest=out/f'boltz_{label}';suffix=f'affinity_{name}.json' if label=='affinity' else f'{name}_model_0.cif'
            script(label,'run_boltz_compatible.py','predict',inp/f'boltz/{name}.yaml','--out_dir',dest,*common,*extra,python=c['boltz_python'],expected=[dest/f'boltz_results_{name}/predictions/{name}/{suffix}'])
        script('analyze','analyze_debian_boltz.py','--source',out,'--output',out/'analysis')
    elif chapter=='15':
        dest=out/'boltzgen_gpu'
        add('configure',[c['boltzgen'],'configure',inp/'boltzgen/penguinpox.yaml','--output',dest,'--protocol','nanobody-anything','--num_designs','2','--budget','2','--steps','design','inverse_folding','folding','analysis','filtering','--num_workers','0','--devices','1','--use_kernels','false','--cache',ROOT/'downloads/huggingface/hub'],expected=[dest/'steps.yaml'])
        add('execute',[c['boltzgen'],'execute',dest],expected=[dest/'final_ranked_designs/final_designs_metrics_2.csv'])
        script('analyze','analyze_debian_boltzgen.py','--source',out,'--output',out/'analysis')
    return tasks

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--chapter',required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
    if not a.config.exists():p.error('Copy config/tools.example.json to config/tools.json and edit your paths first.')
    c=json.loads(a.config.read_text());out=a.output.resolve()
    tasks=plan(a.chapter,c,out)
    if not a.execute:
        print(json.dumps({'execute':False,'chapter':a.chapter,'output':str(out),'tasks':tasks},ensure_ascii=False,indent=2));sys.exit(0)
    if out.exists():raise FileExistsError(f'Use a new output directory: {out}')
    # Check every command before starting any expensive computation.
    for t in tasks:
        executable=t['argv'][0]
        if not (Path(executable).is_file() if '/' in executable else shutil.which(executable)):raise FileNotFoundError(executable)
    out.mkdir(parents=True);(out/'logs').mkdir()
    records=[]
    for t in tasks:
        cwd=Path(t['cwd']);cwd.mkdir(parents=True,exist_ok=True)
        if a.chapter=='12':(out/'runtime/rfantibody').mkdir(parents=True,exist_ok=True)
        env={**os.environ,'MPLBACKEND':'Agg','MPLCONFIGDIR':str(ROOT/'runs/.cache/matplotlib'),'OMP_NUM_THREADS':'4','MKL_NUM_THREADS':'4','XLA_PYTHON_CLIENT_PREALLOCATE':'false',**c.get('extra_env',{}),**t['env']}
        record={**t,'started':datetime.datetime.now().astimezone().isoformat()};start=time.monotonic()
        with (out/'logs'/f"{t['name']}.log").open('w') as log:
            result=subprocess.run(t['argv'],cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
        missing=[x for x in t['expected'] if not Path(x).exists()]
        record.update(returncode=result.returncode,elapsed_seconds=time.monotonic()-start,missing_outputs=missing,execution_passed=result.returncode==0 and not missing)
        records.append(record);(out/'execution.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
        if not record['execution_passed']:raise RuntimeError(f"{t['name']} failed; see {out/'logs'} and execution.json")
        print(f"Completed {t['name']}",flush=True)
    print('Execution completed. Evaluate the scientific metrics and filter flags separately.')
