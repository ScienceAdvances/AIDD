#!/usr/bin/env python3
"""Chapter-oriented entry point. Does not install software or create environments."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
EX=ROOT/'examples/debian13_20261023'
CHAPTERS=json.loads((ROOT/'chapters/index.json').read_text())
OFFLINE=['01','02','03','04','05','06','07','09','10','11','12','13','14','15']

def execute(script,*args):
    command=[sys.executable,str(ROOT/script),*map(str,args)]
    print('+', ' '.join(command),flush=True)
    subprocess.run(command,cwd=ROOT,check=True,env={**os.environ,'MPLBACKEND':'Agg','MPLCONFIGDIR':str(ROOT/'runs/.cache/matplotlib')})

def analyze(chapter,out):
    if out.exists():raise FileExistsError(f'Use a new output directory: {out}')
    if chapter=='01':execute('scripts/analyze_debian13_core.py','--output',out)
    elif chapter=='02':
        shutil.copytree(ROOT/'examples/proteinmpnn',out);execute('scripts/analyze_proteinmpnn.py','--run',out)
    elif chapter in ['03','06','10']:execute('scripts/audit_saved.py',chapter,'--output',out)
    elif chapter=='04':
        execute('scripts/extract_cart_data.py','--data',ROOT/'data/cart','--output',out)
        execute('scripts/plot_cart.py','--input',out/'source_observations.csv','--output',out)
    elif chapter=='05':execute('scripts/build_case_brief.py','--output',out)
    elif chapter=='07':
        out.mkdir(parents=True)
        for mode in ['rescore','none']:
            source=EX/f'gnina_{mode}_gpu';target=out/mode;shutil.copytree(source,target)
            execute('workflows/gnina/gnina_workflow.py','summarize','--run',target,'--output',target/'scores_recomputed.csv')
            for pose in sorted(target.glob('*.sdf')):
                execute('workflows/gnina/gnina_workflow.py','rmsd','--reference',ROOT/'workflows/gnina/example_1iep/reference.sdf','--poses',pose,'--output',target/f'{pose.stem}_rmsd_recomputed.csv')
    elif chapter=='09':execute('scripts/run_rdkit_case.py','--output',out)
    elif chapter=='11':execute('scripts/analyze_debian_rfd.py','--output',out)
    elif chapter=='12':execute('scripts/analyze_debian_rfantibody.py','--output',out)
    elif chapter=='13':execute('scripts/analyze_af3_case.py','--run',EX/'alphafold3_gpu_repaired/human_mdm2_p53_workflow','--output',out)
    elif chapter=='14':execute('scripts/analyze_debian_boltz.py','--output',out)
    elif chapter=='15':execute('scripts/analyze_debian_boltzgen.py','--output',out)
    else:raise ValueError('This is a reading/exercise chapter; see its README.')

def doctor():
    modules=['numpy','pandas','matplotlib','Bio','openpyxl','rdkit','mdtraj','torch']
    tools=['gmx','haddock3','gnina','boltz','boltzgen']
    print(json.dumps({'python':sys.version.split()[0],'modules':{m:importlib.util.find_spec(m) is not None for m in modules},'commands':{m:shutil.which(m) for m in tools}},indent=2))

def verify():
    failures=[];total=0
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        sha,name=line.split('  ',1);p=ROOT/name;total+=1
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=sha:failures.append(name)
    if failures:raise ValueError(f'Changed or missing distributed files: {failures}')
    print(json.dumps({'passed':True,'verified_files':total}))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='action',required=True)
    commands.add_parser('list');commands.add_parser('doctor');commands.add_parser('verify')
    p=commands.add_parser('analyze');p.add_argument('chapter',choices=OFFLINE+['all']);p.add_argument('--output',type=Path,required=True)
    p=commands.add_parser('fetch');p.add_argument('asset',nargs='?',default='list');p.add_argument('--cache-dir',type=Path,default=ROOT/'downloads')
    p=commands.add_parser('run');p.add_argument('chapter',choices=['01','02','03','04','06','07','09','10','11','12','13','14','15']);p.add_argument('--config',type=Path,default=ROOT/'config/tools.json');p.add_argument('--output',type=Path,required=True);p.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.action=='list':
        for c in CHAPTERS:print(c['number'],c['title'],f"chapters/{c['directory']}/README.md")
    elif args.action=='doctor':doctor()
    elif args.action=='verify':verify()
    elif args.action=='fetch':execute('scripts/fetch_assets.py',args.asset,'--cache-dir',args.cache_dir.resolve())
    elif args.action=='run':
        argv=['--chapter',args.chapter,'--config',args.config.resolve(),'--output',args.output.resolve()]
        if args.execute:argv.append('--execute')
        execute('scripts/run_chapter.py',*argv)
    else:
        output=args.output.resolve()
        if args.chapter=='all':
            output.mkdir(parents=True,exist_ok=False)
            for chapter in OFFLINE:analyze(chapter,output/chapter)
        else:analyze(args.chapter,output)
