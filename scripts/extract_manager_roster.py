"""Read private cohort sources and emit a manager roster as JSON.

This script reads existing files only. Its output includes original IDs and is
intended for the private workbook builder, not for publication.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import openpyxl


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "data" / "current" / "metadata"
GENOTYPE = ROOT / "data" / "current" / "genotype"
MAPPING = ROOT / "data" / "processed" / "private" / "dog_id_mapping.tsv"


def normalized_id(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def main() -> None:
    with MAPPING.open(encoding="utf-8", newline="") as handle:
        mapping = {
            row["original_iid"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }

    work = openpyxl.load_workbook(METADATA / "Dataset_Work.xlsx", read_only=True, data_only=True)
    brands: dict[str, set[str]] = defaultdict(set)
    for row in work["Sheet4"].iter_rows(min_row=2, values_only=True):
        if len(row) > 2:
            numeric_id = normalized_id(row[2])
            brand = normalized_id(row[1])
            if numeric_id and brand:
                brands[numeric_id].add(brand)
    work.close()

    template = openpyxl.load_workbook(
        METADATA / "Data Capture Template CTC Canine with DOT Jan 2024.xlsx",
        read_only=True,
        data_only=True,
    )
    outcomes: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in template["Raw Data"].iter_rows(min_row=2, values_only=True):
        if len(row) > 19:
            brand = normalized_id(row[1])
            if brand:
                outcomes[brand].append(
                    {
                        "status": normalized_id(row[3]),
                        "reason": normalized_id(row[19]),
                    }
                )
    template.close()

    records: list[dict[str, str]] = []
    for breed, filename in (("LR", "Battelle_LABs.fam"), ("GSP", "Battelle_GSPs.fam")):
        for line in (GENOTYPE / filename).read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(f"Unexpected FAM record in {filename}")
            original_id = fields[1]
            study_id = mapping[original_id]["study_id"]
            sex = {"1": "Male", "2": "Female", "0": "Unknown"}.get(fields[4], "Unknown")
            phenotype = fields[5]
            outcome = {"1": "Not behaviorally eliminated", "2": "Behaviorally eliminated", "-9": "Unknown"}.get(phenotype)
            if outcome is None:
                raise ValueError(f"Unexpected phenotype {phenotype} for {original_id}")
            linked_brands = sorted(brands.get(original_id, set()))
            if len(linked_brands) > 1:
                raise ValueError(f"Conflicting brand IDs for {original_id}")
            brand = linked_brands[0] if linked_brands else ""
            template_rows = outcomes.get(brand, []) if brand else []
            if len(template_rows) > 1:
                raise ValueError(f"Conflicting template rows for {brand}")
            template_row = template_rows[0] if template_rows else {}
            missing_outcome = phenotype == "-9"
            missing_link = not brand
            need = []
            if missing_outcome:
                need.append("Updated behavioral outcome")
            if missing_link:
                need.append("Brand ID / genotype ID link")
            records.append(
                {
                    "study_id": study_id,
                    "genotype_id": original_id,
                    "brand_id": brand,
                    "breed": breed,
                    "sex": sex,
                    "plink_code": phenotype,
                    "current_outcome": outcome,
                    "template_status": template_row.get("status", ""),
                    "template_elimination_reason": template_row.get("reason", ""),
                    "outcome_missing": "Yes" if missing_outcome else "No",
                    "id_link_missing": "Yes" if missing_link else "No",
                    "information_requested": "; ".join(need),
                }
            )

    if len(records) != 222 or len({row["study_id"] for row in records}) != 222:
        raise ValueError("The roster must contain exactly 222 unique LR/GSP dogs")
    records.sort(key=lambda row: (row["outcome_missing"] != "Yes", row["id_link_missing"] != "Yes", row["breed"], row["study_id"]))
    print(json.dumps(records, separators=(",", ":")))


if __name__ == "__main__":
    main()
