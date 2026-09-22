"""Create pseudonymous PLINK datasets for the Labrador and GSP cohorts.

The source PLINK files are treated as read-only. Genotype and variant files are
copied byte-for-byte; only sample identifiers in the FAM files are replaced.
The private crosswalk is written beneath data/processed/private, which is
excluded from version control along with all study data.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DogId:
    original_fid: str
    original_iid: str
    study_id: str
    breed: str


COHORTS = {
    "LR": "Battelle_LABs",
    "GSP": "Battelle_GSPs",
}
COMBINED_SOURCE = "Battelle_GSPLABs"


def read_fam(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if any(len(row) != 6 for row in rows):
        raise ValueError(f"Expected six columns in {path}")
    return rows


def build_mapping(genotype_dir: Path) -> list[DogId]:
    mapping: list[DogId] = []
    seen: set[str] = set()

    for breed, stem in COHORTS.items():
        rows = read_fam(genotype_dir / f"{stem}.fam")
        width = max(3, len(str(len(rows))))
        for index, row in enumerate(rows, start=1):
            original_fid, original_iid = row[0], row[1]
            if original_iid in seen:
                raise ValueError(f"Dog {original_iid} occurs in more than one breed cohort")
            seen.add(original_iid)
            mapping.append(
                DogId(
                    original_fid=original_fid,
                    original_iid=original_iid,
                    study_id=f"{breed}{index:0{width}d}",
                    breed=breed,
                )
            )

    return mapping


def write_mapping(path: Path, mapping: list[DogId]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "breed", "original_fid", "original_iid"])
        for dog in mapping:
            writer.writerow([dog.study_id, dog.breed, dog.original_fid, dog.original_iid])


def write_pseudonymous_dataset(
    source_dir: Path,
    output_dir: Path,
    source_stem: str,
    output_stem: str,
    by_iid: dict[str, DogId],
) -> None:
    source_fam = read_fam(source_dir / f"{source_stem}.fam")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_fam = output_dir / f"{output_stem}.fam"
    with output_fam.open("w", encoding="utf-8", newline="") as handle:
        for row in source_fam:
            dog = by_iid.get(row[1])
            if dog is None:
                raise ValueError(f"No pseudonymous ID for {row[1]} in {source_stem}")
            row[0] = dog.study_id
            row[1] = dog.study_id
            handle.write(" ".join(row) + "\n")

    for extension in ("bed", "bim"):
        shutil.copyfile(
            source_dir / f"{source_stem}.{extension}",
            output_dir / f"{output_stem}.{extension}",
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    source_dir = project_root / "data" / "current" / "genotype"
    processed_dir = project_root / "data" / "processed"
    output_dir = processed_dir / "genotype"

    mapping = build_mapping(source_dir)
    by_iid = {dog.original_iid: dog for dog in mapping}
    write_mapping(processed_dir / "private" / "dog_id_mapping.tsv", mapping)

    for breed, source_stem in COHORTS.items():
        write_pseudonymous_dataset(
            source_dir,
            output_dir,
            source_stem,
            f"2022TSA_{breed}",
            by_iid,
        )

    write_pseudonymous_dataset(
        source_dir,
        output_dir,
        COMBINED_SOURCE,
        "2022TSA_LR_GSP",
        by_iid,
    )

    print(f"Created mapping for {len(mapping)} dogs")
    print(f"Private mapping: {processed_dir / 'private' / 'dog_id_mapping.tsv'}")
    print(f"Pseudonymous PLINK datasets: {output_dir}")


if __name__ == "__main__":
    main()
