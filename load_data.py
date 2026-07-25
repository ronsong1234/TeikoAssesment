"""Create and populate the SQLite database for the cell-count dataset.

Run from any working directory with:

    python load_data.py

The input CSV and output database are resolved relative to this script so the
loader behaves consistently regardless of the caller's current directory.
"""

from __future__ import annotations

import csv
import os
import sqlite3
from pathlib import Path
from typing import Iterator


ROOT_DIR = Path(__file__).resolve().parent
CSV_PATH = ROOT_DIR / "cell-count.csv"
DATABASE_PATH = ROOT_DIR / "cell-count.db"
TEMP_DATABASE_PATH = ROOT_DIR / "cell-count.db.tmp"

CELL_COLUMNS = (
    ("b_cell", "B cell"),
    ("cd8_t_cell", "CD8 T cell"),
    ("cd4_t_cell", "CD4 T cell"),
    ("nk_cell", "NK cell"),
    ("monocyte", "Monocyte"),
)

EXPECTED_COLUMNS = (
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *(column for column, _ in CELL_COLUMNS),
)

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
        CHECK (length(trim(project_id)) > 0)
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY
        CHECK (length(trim(subject_id)) > 0),
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition TEXT NOT NULL CHECK (length(trim(condition)) > 0),
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL CHECK (sex IN ('F', 'M')),
    treatment TEXT NOT NULL CHECK (length(trim(treatment)) > 0),
    response TEXT CHECK (response IN ('yes', 'no') OR response IS NULL)
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY
        CHECK (length(trim(sample_id)) > 0),
    subject_id TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type TEXT NOT NULL CHECK (length(trim(sample_type)) > 0),
    time_from_treatment_start INTEGER NOT NULL
        CHECK (time_from_treatment_start >= 0),
    UNIQUE (subject_id, sample_type, time_from_treatment_start)
);

CREATE TABLE cell_types (
    cell_type_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    display_name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL REFERENCES samples(sample_id),
    cell_type_id INTEGER NOT NULL REFERENCES cell_types(cell_type_id),
    cell_count INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, cell_type_id)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_subjects_analysis
    ON subjects(condition, treatment, response, sex);
CREATE INDEX idx_samples_subject_time
    ON samples(subject_id, time_from_treatment_start);
CREATE INDEX idx_cell_counts_type ON cell_counts(cell_type_id);

CREATE VIEW sample_cell_counts AS
SELECT
    sub.project_id AS project,
    sub.subject_id AS subject,
    sub.condition,
    sub.age,
    sub.sex,
    sub.treatment,
    sub.response,
    sam.sample_id AS sample,
    sam.sample_type,
    sam.time_from_treatment_start,
    MAX(CASE WHEN ct.name = 'b_cell' THEN cc.cell_count END) AS b_cell,
    MAX(CASE WHEN ct.name = 'cd8_t_cell' THEN cc.cell_count END) AS cd8_t_cell,
    MAX(CASE WHEN ct.name = 'cd4_t_cell' THEN cc.cell_count END) AS cd4_t_cell,
    MAX(CASE WHEN ct.name = 'nk_cell' THEN cc.cell_count END) AS nk_cell,
    MAX(CASE WHEN ct.name = 'monocyte' THEN cc.cell_count END) AS monocyte
FROM samples AS sam
JOIN subjects AS sub ON sub.subject_id = sam.subject_id
JOIN cell_counts AS cc ON cc.sample_id = sam.sample_id
JOIN cell_types AS ct ON ct.cell_type_id = cc.cell_type_id
GROUP BY sam.sample_id;

CREATE VIEW cell_population_frequencies AS
SELECT
    cc.sample_id AS sample,
    SUM(cc.cell_count) OVER (PARTITION BY cc.sample_id) AS total_count,
    ct.name AS population,
    cc.cell_count AS count,
    100.0 * cc.cell_count
        / NULLIF(
            SUM(cc.cell_count) OVER (PARTITION BY cc.sample_id),
            0
        ) AS percentage
FROM cell_counts AS cc
JOIN cell_types AS ct ON ct.cell_type_id = cc.cell_type_id;
"""


def required_text(row: dict[str, str], column: str, row_number: int) -> str:
    """Return a non-empty, trimmed field value."""
    value = (row.get(column) or "").strip()
    if not value:
        raise ValueError(f"Row {row_number}: {column!r} must not be empty")
    return value


def nonnegative_integer(
    row: dict[str, str], column: str, row_number: int
) -> int:
    """Parse a field as a non-negative integer with a useful error message."""
    value = required_text(row, column, row_number)
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(
            f"Row {row_number}: {column!r} must be an integer, got {value!r}"
        ) from exc
    if parsed < 0:
        raise ValueError(
            f"Row {row_number}: {column!r} must be non-negative, got {parsed}"
        )
    return parsed


def rows_from_csv() -> Iterator[tuple[int, dict[str, str]]]:
    """Yield validated CSV records and their one-based file row numbers."""
    if not CSV_PATH.is_file():
        raise FileNotFoundError(
            f"Input file not found: {CSV_PATH}\n"
            "Place cell-count.csv beside load_data.py and run the script again."
        )

    with CSV_PATH.open("r", encoding="utf-8-sig", newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        actual_columns = tuple(reader.fieldnames or ())
        if actual_columns != EXPECTED_COLUMNS:
            raise ValueError(
                "Unexpected CSV columns.\n"
                f"Expected: {', '.join(EXPECTED_COLUMNS)}\n"
                f"Found:    {', '.join(actual_columns)}"
            )
        for row_number, row in enumerate(reader, start=2):
            yield row_number, row


def load_database(connection: sqlite3.Connection) -> tuple[int, int, int]:
    """Create the schema and insert all CSV records in one transaction."""
    connection.executescript(SCHEMA)
    connection.executemany(
        "INSERT INTO cell_types (cell_type_id, name, display_name) VALUES (?, ?, ?)",
        (
            (cell_type_id, column, display_name)
            for cell_type_id, (column, display_name) in enumerate(
                CELL_COLUMNS, start=1
            )
        ),
    )

    projects: set[str] = set()
    subjects: dict[str, tuple[object, ...]] = {}
    sample_count = 0

    for row_number, row in rows_from_csv():
        project = required_text(row, "project", row_number)
        subject = required_text(row, "subject", row_number)
        response_value = (row.get("response") or "").strip().lower()
        response = response_value or None
        if response not in (None, "yes", "no"):
            raise ValueError(
                f"Row {row_number}: 'response' must be yes, no, or empty"
            )

        subject_record = (
            subject,
            project,
            required_text(row, "condition", row_number),
            nonnegative_integer(row, "age", row_number),
            required_text(row, "sex", row_number).upper(),
            required_text(row, "treatment", row_number),
            response,
        )
        previous_record = subjects.get(subject)
        if previous_record is not None and previous_record != subject_record:
            raise ValueError(
                f"Row {row_number}: inconsistent metadata for subject {subject!r}"
            )

        if project not in projects:
            connection.execute(
                "INSERT INTO projects (project_id) VALUES (?)", (project,)
            )
            projects.add(project)

        if previous_record is None:
            connection.execute(
                """
                INSERT INTO subjects (
                    subject_id, project_id, condition, age, sex, treatment, response
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                subject_record,
            )
            subjects[subject] = subject_record

        sample = required_text(row, "sample", row_number)
        connection.execute(
            """
            INSERT INTO samples (
                sample_id, subject_id, sample_type, time_from_treatment_start
            ) VALUES (?, ?, ?, ?)
            """,
            (
                sample,
                subject,
                required_text(row, "sample_type", row_number),
                nonnegative_integer(
                    row, "time_from_treatment_start", row_number
                ),
            ),
        )
        connection.executemany(
            """
            INSERT INTO cell_counts (sample_id, cell_type_id, cell_count)
            VALUES (?, ?, ?)
            """,
            (
                (
                    sample,
                    cell_type_id,
                    nonnegative_integer(row, column, row_number),
                )
                for cell_type_id, (column, _) in enumerate(
                    CELL_COLUMNS, start=1
                )
            ),
        )
        sample_count += 1

    if sample_count == 0:
        raise ValueError("The CSV contains no data rows")

    return len(projects), len(subjects), sample_count


def main() -> None:
    """Build a fresh database and atomically replace any previous output."""
    TEMP_DATABASE_PATH.unlink(missing_ok=True)
    connection: sqlite3.Connection | None = None
    try:
        connection = sqlite3.connect(TEMP_DATABASE_PATH)
        try:
            connection.execute("PRAGMA foreign_keys = ON")
            project_count, subject_count, sample_count = load_database(connection)
            integrity_result = connection.execute("PRAGMA integrity_check").fetchone()
            if integrity_result != ("ok",):
                raise RuntimeError(
                    f"SQLite integrity check failed: {integrity_result}"
                )
            connection.commit()
        finally:
            connection.close()
            connection = None
        os.replace(TEMP_DATABASE_PATH, DATABASE_PATH)
    except Exception:
        if connection is not None:
            connection.close()
        TEMP_DATABASE_PATH.unlink(missing_ok=True)
        raise

    print(
        f"Created {DATABASE_PATH.name}: "
        f"{project_count} projects, {subject_count} subjects, "
        f"{sample_count} samples, {sample_count * len(CELL_COLUMNS)} cell counts."
    )


if __name__ == "__main__":
    main()
