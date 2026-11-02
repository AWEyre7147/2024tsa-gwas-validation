"""Check notebook helper loading in isolated processes, including reordered imports."""
from pathlib import Path
import ast
import json
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


class NotebookImportTests(unittest.TestCase):
    def test_fresh_kernel_and_organized_imports(self):
        for path in sorted((ROOT / 'notebooks').glob('*.ipynb')):
            book = json.loads(path.read_text(encoding='utf-8'))
            code = [c for c in book['cells']
                    if c['cell_type'] == 'code' and ''.join(c['source']).strip()]
            if not code:
                self.assertGreaterEqual(int(path.name[:2]), 11, path.name)
                continue
            cell = code[0]
            source = ''.join(cell['source'])
            body = ast.parse(source).body
            helper_nodes = [n for n in body if isinstance(n, ast.Assign)
                            and isinstance(n.value, ast.Attribute)
                            and isinstance(n.value.value, ast.Call)
                            and ast.unparse(n.value.value.func) == 'importlib.import_module']
            self.assertTrue(helper_nodes, path.name)
            end = max(n.end_lineno for n in helper_nodes)
            startup = ast.parse('\n'.join(source.splitlines()[:end]))
            # Simulate an editor moving ordinary imports ahead of all other code.
            imports = [n for n in startup.body if isinstance(n, (ast.Import, ast.ImportFrom))
                       and not (isinstance(n, ast.ImportFrom) and n.module.startswith('IPython'))]
            other = [n for n in startup.body if not isinstance(n, (ast.Import, ast.ImportFrom))]
            organized = ast.unparse(ast.Module(body=imports + other, type_ignores=[]))
            check = "\nimport notebook_support, notebook_review\nassert callable(require_saved)\n"
            check += f"assert Path(notebook_support.__file__).resolve() == Path({str(ROOT / 'scripts/notebook_support.py')!r}).resolve()\n"
            for folder in (ROOT, ROOT / 'notebooks'):
                with self.subTest(notebook=path.name, folder=folder.name):
                    result = subprocess.run([sys.executable, '-I', '-c', organized + check],
                                            cwd=folder, text=True, capture_output=True)
                    self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
