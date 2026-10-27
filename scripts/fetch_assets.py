"""Download one explicitly selected, checksum-locked teaching asset."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request
import zipfile
ROOT=Path(__file__).resolve().parents[1]
ASSETS=json.loads((ROOT/'config/assets.json').read_text())
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def fetch(name,cache):
    if name not in ASSETS:raise ValueError(f'Unknown asset: {name}; use list')
    a=ASSETS[name];cache.mkdir(parents=True,exist_ok=True);cached=cache/(name+('.zip' if a['kind']=='zip' else '.bin'))
    if not cached.exists():
        partial=cached.with_suffix(cached.suffix+'.partial')
        req=urllib.request.Request(a['url'],headers={'User-Agent':'AIDD-practice/1.0'})
        with urllib.request.urlopen(req,timeout=120) as src,partial.open('wb') as dst:shutil.copyfileobj(src,dst)
        if partial.stat().st_size!=a['bytes'] or digest(partial)!=a['sha256']:
            raise ValueError('Downloaded asset differs from the Debian snapshot. Kept .partial for inspection; did not install it.')
        partial.replace(cached)
    if cached.stat().st_size!=a['bytes'] or digest(cached)!=a['sha256']:raise ValueError('Cached asset is damaged or changed')
    dest=ROOT/a['destination']
    if a['kind']=='file':
        dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists() and digest(dest)!=a['sha256']:raise FileExistsError(f'Different existing file: {dest}')
        shutil.copyfile(cached,dest)
    else:
        dest.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(cached) as archive:
            members=[]
            for item in archive.infolist():
                path=Path(item.filename)
                if path.is_absolute() or '..' in path.parts or ((item.external_attr>>16)&0o170000)==0o120000:raise ValueError('Unsafe archive member')
                target=(dest/path).resolve()
                if not target.is_relative_to(dest.resolve()):raise ValueError('Archive escapes destination')
                if item.is_dir():continue
                # Compare existing content before changing any input files.
                if target.exists() and hashlib.sha256(archive.read(item)).hexdigest()!=digest(target):raise FileExistsError(f'Different input already exists: {target}')
                members.append((item,target))
            for item,target in members:
                target.parent.mkdir(parents=True,exist_ok=True)
                with archive.open(item) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
    print(json.dumps(dict(asset=name,sha256=a['sha256'],destination=str(dest),verified=True),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('asset',nargs='?',default='list');p.add_argument('--cache-dir',type=Path,default=ROOT/'downloads');args=p.parse_args()
    if args.asset=='list':
        for name,a in ASSETS.items():print(name,a['description'])
    else:fetch(args.asset,args.cache_dir.resolve())
