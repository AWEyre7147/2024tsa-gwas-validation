"""Check the Git publication boundary, notebook sources and local documentation links."""
from pathlib import Path
import ast,json,re,subprocess,sys
from urllib.parse import unquote
from build_notebook_previews import check_outputs, preview_record

ROOT=Path(__file__).resolve().parents[1]
PUBLIC_DIRS={'notebooks','scripts','config','docs'}
ROOT_FILES={'.gitignore','.gitattributes','README.md'}
PRIVATE_PARTS={'data','results','.venv','.virtual_documents','.ipynb_checkpoints','__pycache__','.obsidian','raw','node_modules'}
SECRET_PATTERNS=[re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
                 re.compile(r'\b(?:ghp_|github_pat_|sk-proj-)[A-Za-z0-9_]{20,}'),
                 re.compile(r'\bAKIA[A-Z0-9]{16}\b')]

def git_paths(staged):
    args=['git','diff','--cached','--name-only','--diff-filter=ACMR','-z'] if staged else ['git','ls-files','--cached','--others','--exclude-standard','-z']
    return sorted(set(subprocess.check_output(args,cwd=ROOT).decode('utf-8').split('\0'))-{''})

def check_notebook_contents(staged=False):
    """Require the overview to match chapter titles and Markdown headings in order."""
    sections=['## Table of Contents\n'];actual=None
    try:
        names=sorted(name for name in git_paths(staged=False)
                     if name.startswith('notebooks/') and name.endswith('.ipynb')
                     and (ROOT/name).exists())
        if staged:
            names=sorted(subprocess.check_output(
                ['git','ls-files','notebooks/*.ipynb'],cwd=ROOT).decode().splitlines())
        for name in names:
            text=(subprocess.check_output(['git','show',':'+name],cwd=ROOT).decode('utf-8')
                  if staged else (ROOT/name).read_text(encoding='utf-8'))
            book=json.loads(text);title=''.join(book['cells'][0]['source']).strip().removeprefix('# ')
            sections.append(f'### {Path(name).stem} — {title}\n')
            headings=[]
            for cell in book['cells'][1:]:
                if cell['cell_type']!='markdown':continue
                source=''.join(cell['source'])
                if source.startswith('## Table of Contents'):
                    actual=source;headings.append('Table of Contents');continue
                headings.extend(re.findall(r'^#{2,6} (.+?)\s*$',source,re.M))
            if headings:sections.append('\n'.join('- '+heading for heading in headings)+'\n')
        if actual is None or actual.strip()!='\n'.join(sections).strip():
            return ['Notebook 00 Table of Contents does not match current chapter titles and section headings.']
    except (OSError,ValueError,KeyError,IndexError,subprocess.CalledProcessError):
        return ['Unable to verify notebook Table of Contents.']
    return []

def check(staged=False):
    errors=[];books=0;code_cells=0;public=[];notes={p.stem for p in (ROOT/'docs').rglob('*.md') if 'raw' not in p.parts}
    preview_path='config/notebook_previews.json'
    try:
        preview_text=(subprocess.check_output(['git','show',':'+preview_path],cwd=ROOT).decode('utf-8')
                      if staged else (ROOT/preview_path).read_text(encoding='utf-8'))
        previews=json.loads(preview_text)['notebooks']
    except (OSError,ValueError,KeyError,subprocess.CalledProcessError):
        previews={}
    for name in git_paths(staged):
        relative=Path(name);path=ROOT/relative
        if not path.exists():continue
        if relative.parts[0] not in PUBLIC_DIRS and name not in ROOT_FILES:errors.append(f'Unexpected public path: {name}');continue
        if any(part in PRIVATE_PARTS for part in relative.parts) or '.private.' in name or relative.suffix.lower() in {'.bed','.bim','.fam','.vcf','.bcf','.rds','.gz','.zip','.exe','.png','.svg','.pdf'}:
            errors.append(f'Private/generated artifact in Git selection: {name}');continue
        if not path.resolve().is_relative_to(ROOT):errors.append(f'External symlink target: {name}');continue
        if path.stat().st_size>10_000_000:errors.append(f'Unexpectedly large public source: {name}');continue
        try:
            text=(subprocess.check_output(['git','show',':'+name],cwd=ROOT).decode('utf-8-sig')
                  if staged else path.read_text(encoding='utf-8-sig'))
        except UnicodeDecodeError:errors.append(f'Non-text public artifact: {name}');continue
        if any(pattern.search(text) for pattern in SECRET_PATTERNS):errors.append(f'Potential credential material: {name}')
        public.append(name)
        if path.suffix=='.py':
            try:ast.parse(text,filename=name)
            except SyntaxError as exc:errors.append(f'Python syntax: {name}:{exc.lineno}')
        if path.suffix=='.json':
            try:json.loads(text)
            except ValueError:errors.append(f'Invalid JSON: {name}')
        documents=[]
        if path.suffix=='.md':documents=[text]
        if path.suffix=='.ipynb':
            books+=1
            try:book=json.loads(text)
            except ValueError:errors.append(f'Invalid notebook JSON: {name}');continue
            if book.get('nbformat')!=4:errors.append(f'Unexpected notebook format: {name}')
            if book.get('metadata',{}).get('widgets'):errors.append(f'Embedded widget state: {name}')
            has_outputs=any(c.get('outputs') or c.get('execution_count') is not None for c in book.get('cells',[]))
            if has_outputs:
                if previews.get(path.name)!=preview_record(book):
                    errors.append(f'Notebook preview differs from reviewed output/source hashes: {name}')
                try:check_outputs(book)
                except ValueError as exc:errors.append(f'Unsafe notebook output: {name}: {exc}')
            for i,cell in enumerate(book.get('cells',[]),1):
                source=''.join(cell['source'])
                if cell.get('attachments'):errors.append(f'Embedded notebook attachment: {name}:{i}')
                if cell['cell_type']=='markdown':documents.append(source)
                if cell['cell_type']=='code':
                    code_cells+=1
                    try:compile(source,name,'exec')
                    except SyntaxError as exc:errors.append(f'Notebook code syntax: {name}:{i}:{exc.lineno}')
        for document in documents:
            for target in re.findall(r'\[[^\]\n]*\]\(([^\n]*?)\)',document):
                target=target.strip().strip('<>')
                if target.startswith(('http:','https:','mailto:','#')):continue
                if target.startswith(('attachment:','data:')):errors.append(f'Embedded/missing image reference: {name}');continue
                resolved=path.parent/unquote(target.split('#',1)[0])
                if not resolved.exists():errors.append(f'Broken local link: {name} -> {target}')
            for target in re.findall(r'\[\[([^\]]+)\]\]',document):
                stem=target.split('|',1)[0].split('#',1)[0].removesuffix('.md')
                if stem and stem not in notes:errors.append(f'Broken Obsidian link: {name} -> {stem}')
    if not staged and books!=13:errors.append(f'Expected 13 canonical notebooks, found {books}')
    errors.extend(check_notebook_contents(staged))
    for private in ['data/01-input/manifest.private.json','results/private/version2/notebooks/03_structure_and_model_design.ipynb','.virtual_documents/example.py']:
        if subprocess.run(['git','check-ignore','-q',private],cwd=ROOT).returncode:errors.append(f'Private path is not ignored: {private}')
    if errors:
        print('\n'.join(errors));return 1
    print(f'Passed: {len(public)} public sources, {books} notebooks, {code_cells} code cells; reviewed output hashes, static output formats, local links, syntax and private-path exclusions checked.')
    print('This check inspects current files; it does not audit all historical Git objects or replace review of intended public content.')
    return 0

if __name__=='__main__':sys.exit(check('--staged' in sys.argv))
