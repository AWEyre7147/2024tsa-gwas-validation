"""Execute saved-result review and embed static outputs for GitHub (no new analyses).

Run in the study's scientific Python environment. Local authorized inputs and
completed result files are required. This command does not publish to GitHub.
"""
from pathlib import Path
from notebook_review import check_review_code
import argparse
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = ROOT / "results/private/notebook-readability-20261025"
ALLOWED_MIME = {"text/plain", "text/html", "text/markdown", "image/png"}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def source_digest(book):
    return digest([{k: c.get(k) for k in ("cell_type", "id", "source")} for c in book["cells"]])


def preview_record(book):
    return {"source_sha256": source_digest(book), "cells": {
        c["id"]: {"execution_count": c.get("execution_count"),
                  "outputs_sha256": digest(c.get("outputs", []))}
        for c in book["cells"] if c["cell_type"] == "code"
    }}


def check_outputs(book, identifiers=()):
    """Reject active content, errors, local user paths and known dog identifiers."""
    for cell in book["cells"]:
        for output in cell.get("outputs", []):
            if output["output_type"] not in {"stream", "display_data", "execute_result"}:
                raise ValueError("Only static successful outputs may appear in the preview")
            data = output.get("data", {})
            if set(data) - ALLOWED_MIME:
                raise ValueError(f"Unsupported output formats: {set(data) - ALLOWED_MIME}")
            text = json.dumps({k: v for k, v in data.items() if k != "image/png"}, ensure_ascii=False)
            text += json.dumps(output.get("text", ""), ensure_ascii=False)
            if re.search(r"(?:/home/|/mnt/[a-z]/Users/|[A-Z]:[\\/]+Users[\\/])", text, re.I):
                raise ValueError("A local user path appears in a notebook output")
            if re.search(r"<\s*(?:script|iframe|object|embed)\b|javascript:|on(?:load|error)\s*=", text, re.I):
                raise ValueError("Active HTML is not allowed in notebook outputs")
            if any(re.search(r"(?<![\w])" + re.escape(identifier) + r"(?![\w])", text)
                   for identifier in identifiers):
                raise ValueError("A private dog identifier appears in notebook output")


def known_identifiers():
    identifiers = set()
    for path in (ROOT / "data/01-input").rglob("*.fam"):
        for line in path.read_text().splitlines():
            columns = line.split()
            if len(columns) >= 2 and len(columns[1]) >= 5:
                identifiers.add(columns[1])
    return identifiers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", required=True,
                        help="Execute saved-result review without launching analysis commands")
    args = parser.parse_args()
    import nbformat
    from nbclient import NotebookClient
    from nbconvert import HTMLExporter
    from jupyter_client import KernelManager
    from jupyter_client.kernelspec import KernelSpec

    PRIVATE.mkdir(parents=True, exist_ok=True)
    identifiers = known_identifiers()
    manifest = {"description": "Reviewed static study outputs; new analyses disabled",
                "notebooks": {}}
    for path in sorted((ROOT / "notebooks").glob("*.ipynb")):
        book = nbformat.read(path, as_version=4)
        check_review_code(book)
        for cell in book.cells:
            if cell.cell_type != "code":
                continue
            cell.outputs = []
            cell.execution_count = None
        print("Reviewing saved results:", path.name, flush=True)
        manager = KernelManager(kernel_name="python3")
        manager._kernel_spec = KernelSpec(argv=[sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
                                         display_name="Study Python", language="python")
        client = NotebookClient(book, timeout=1800, kernel_name="python3", km=manager,
                                resources={"metadata": {"path": str(ROOT / "notebooks")}})
        try:
            client.execute()
        finally:
            if manager.has_kernel:
                manager.shutdown_kernel(now=True)
        nbformat.write(book, PRIVATE / path.name)
        check_outputs(book, identifiers)
        book.metadata["study_preview"] = {"mode": "saved-result review",
                                          "new_analyses_enabled": False}
        # Match nbformat's list-of-lines representation when computing receipts.
        serialized = nbformat.writes(book)
        public_book = json.loads(serialized)
        manifest["notebooks"][path.name] = preview_record(public_book)
        path.write_text(serialized, encoding="utf-8")
        exporter = HTMLExporter(template_name="classic")
        exporter.exclude_input_prompt = True
        html, _ = exporter.from_notebook_node(book)
        (PRIVATE / (path.stem + ".html")).write_text(html, encoding="utf-8")
    (ROOT / "config/notebook_previews.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("Saved 13 notebook previews and private HTML review copies. No new analyses were enabled.")


if __name__ == "__main__":
    main()
