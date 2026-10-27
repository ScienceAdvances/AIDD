"""Validate the portable public payload and optionally regenerate SHA256SUMS."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote
ROOT=Path(__file__).resolve().parents[1]
SKIP={'runs','downloads','optional_data','.git','.pytest_cache','__pycache__'}
WEIGHTS={'.pt','.pth','.ckpt','.safetensors','.sif','.bin','.pyc'}
TEXT={'.md','.py','.json','.jsonl','.yaml','.yml','.cfg','.tsv','.csv','.pdb','.cif','.sdf','.txt','.fa','.fasta','.mdp','.toml'}
def payload():
    for p in sorted(ROOT.rglob('*')):
        rel=p.relative_to(ROOT)
        if any(x in SKIP or x.startswith('.venv') for x in rel.parts):continue
        if rel.as_posix()=='config/tools.json' or p.suffix in WEIGHTS or p.name.startswith('.env'):continue
        if rel.as_posix().startswith('data/md/md-intro-tutorial-main/data/reference/'):continue
        if p.is_symlink():raise ValueError(f'No symlinks in the distribution: {rel}')
        if p.is_file():yield p

def check():
    files=list(payload());badlinks=[];pythonfiles=0
    for p in files:
        rel=p.relative_to(ROOT)
        if p.stat().st_size>=25*2**20:raise ValueError(f'File too large for browser upload: {rel}')
        if p.suffix in TEXT:
            s=p.read_text()
            if re.search('/home/'+'tim/|/User'+'s/',s):raise ValueError(f'Author-specific path: {rel}')
        if p.suffix=='.py':ast.parse(p.read_text(),filename=str(rel));pythonfiles+=1
        if p.suffix=='.md':
            s=re.sub(r'(?ms)^```[^\n]*\n.*?^```[^\n]*','',p.read_text())
            for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',s):
                target=target.strip('<>');target=unquote(target.split('#')[0])
                if not target or re.match(r'^[a-zA-Z][\w+.-]*:',target):continue
                dest=p.parent/target
                if not dest.exists():badlinks.append((str(rel),target))
    if badlinks:raise ValueError('Missing local Markdown targets:\n'+'\n'.join(f'{p}: {t}' for p,t in badlinks))
    return dict(passed=True,files=len(files),bytes=sum(p.stat().st_size for p in files),python_syntax_files=pythonfiles,local_markdown_links='passed',author_absolute_paths='none',symlinks='none',max_file_bytes=max(p.stat().st_size for p in files))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-manifest',action='store_true');a=p.parse_args()
    result=check()
    if a.write_manifest:
        entries=[]
        for f in payload():
            if f.name=='SHA256SUMS' and f.parent==ROOT:continue
            entries.append(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(ROOT).as_posix()}")
        (ROOT/'SHA256SUMS').write_text('\n'.join(entries)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
