"""Checks for notebooks that read saved results and show analysis commands as examples."""
from pathlib import Path
import ast


def require_saved(path):
    """Explain missing saved results before a notebook attempts new computation."""
    if not Path(path).exists():
        raise FileNotFoundError(
            f"Saved results are required for review: {Path(path).name}. "
            "The published notebook already includes its recorded results. "
            "New calculations use the corresponding analysis scripts with "
            "the required inputs, tools and output folders."
        )


def check_review_code(book):
    """Reject executable external analysis commands in the saved-result notebooks."""
    for cell in book['cells']:
        if cell['cell_type'] != 'code':
            continue
        for node in ast.walk(ast.parse(''.join(cell['source']))):
            if isinstance(node, ast.Call) and ast.unparse(node.func) in ('run_script', 'subprocess.run'):
                raise ValueError('Analysis launch commands must remain commented in the review notebooks')
