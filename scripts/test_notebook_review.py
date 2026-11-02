"""Check that saved-result notebooks cannot launch external analyses accidentally."""
from pathlib import Path
import json
import unittest
from unittest.mock import patch

from notebook_review import require_saved, check_review_code


class NotebookReviewTests(unittest.TestCase):
    def test_missing_results_stop_review(self):
        with patch('notebook_review.Path.exists', return_value=False):
            with self.assertRaisesRegex(FileNotFoundError, 'Saved results'):
                require_saved(Path('missing-review/manifest.private.json'))

    def test_active_analysis_commands_are_rejected(self):
        for source in ['run_script(root, "imputation_workflow.py")', 'subprocess.run(command)']:
            with self.subTest(source=source), self.assertRaises(ValueError):
                check_review_code({'cells': [{'cell_type': 'code', 'source': source}]})

    def test_commented_commands_are_allowed(self):
        check_review_code({'cells': [{'cell_type': 'code', 'source': '# run_script(root, "example.py")\nvalue = 1'}]})

    def test_canonical_notebooks_only_review_saved_results(self):
        root = Path(__file__).resolve().parents[1]
        for path in (root/'notebooks').glob('*.ipynb'):
            with self.subTest(notebook=path.name):
                check_review_code(json.loads(path.read_text(encoding='utf-8')))


if __name__ == '__main__':
    unittest.main()
