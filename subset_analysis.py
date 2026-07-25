"""Analyze baseline melanoma PBMC samples from miraclib-treated patients."""

import csv
import sqlite3
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent
DATABASE_PATH = ROOT_DIR / "cell-count.db"
SUBSET_PATH = ROOT_DIR / "baseline_melanoma_pbmc_miraclib.csv"

BASE_FILTER = """
FROM samples AS sample
JOIN subjects AS subject ON subject.subject_id = sample.subject_id
WHERE subject.condition = 'melanoma'
  AND subject.treatment = 'miraclib'
  AND sample.sample_type = 'PBMC'
  AND sample.time_from_treatment_start = 0
"""

SAMPLE_QUERY = f"""
SELECT
    subject.project_id AS project,
    subject.subject_id AS subject,
    subject.response,
    subject.sex,
    sample.sample_id AS sample,
    sample.sample_type,
    sample.time_from_treatment_start
{BASE_FILTER}
ORDER BY subject.project_id, subject.subject_id, sample.sample_id
"""

PROJECT_QUERY = f"""
SELECT subject.project_id AS project, COUNT(*) AS sample_count
{BASE_FILTER}
GROUP BY subject.project_id
ORDER BY subject.project_id
"""

RESPONSE_QUERY = f"""
SELECT subject.response, COUNT(DISTINCT subject.subject_id) AS subject_count
{BASE_FILTER}
GROUP BY subject.response
ORDER BY subject.response
"""

SEX_QUERY = f"""
SELECT subject.sex, COUNT(DISTINCT subject.subject_id) AS subject_count
{BASE_FILTER}
GROUP BY subject.sex
ORDER BY subject.sex
"""


def query_database():
    if not DATABASE_PATH.exists():
        raise FileNotFoundError(
            "cell-count.db was not found. Run 'python load_data.py' first."
        )

    connection = sqlite3.connect(DATABASE_PATH)
    connection.row_factory = sqlite3.Row
    try:
        samples = connection.execute(SAMPLE_QUERY).fetchall()
        projects = connection.execute(PROJECT_QUERY).fetchall()
        responses = connection.execute(RESPONSE_QUERY).fetchall()
        sexes = connection.execute(SEX_QUERY).fetchall()
    finally:
        connection.close()

    return samples, projects, responses, sexes


def write_subset(samples):
    columns = (
        "project",
        "subject",
        "response",
        "sex",
        "sample",
        "sample_type",
        "time_from_treatment_start",
    )
    with SUBSET_PATH.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(columns)
        writer.writerows(tuple(row[column] for column in columns) for row in samples)


def display_samples(samples):
    print("Qualifying baseline samples")
    print(
        f"{'project':<8} {'subject':<10} {'response':<10} "
        f"{'sex':<5} {'sample':<14}"
    )
    print("-" * 51)
    for row in samples:
        print(
            f"{row['project']:<8} {row['subject']:<10} "
            f"{row['response']:<10} {row['sex']:<5} {row['sample']:<14}"
        )


def display_summary(title, rows, label_column, count_column):
    print(f"\n{title}")
    print(f"{'group':<16} {'count':>8}")
    print("-" * 25)
    for row in rows:
        print(f"{row[label_column]:<16} {row[count_column]:>8}")


def main():
    samples, projects, responses, sexes = query_database()
    write_subset(samples)
    display_samples(samples[:20])
    print(f"\nShowing 20 of {len(samples)} qualifying samples.")
    display_summary("Samples by project", projects, "project", "sample_count")
    display_summary(
        "Subjects by response", responses, "response", "subject_count"
    )
    display_summary("Subjects by sex", sexes, "sex", "subject_count")
    print(f"\nTotal qualifying samples: {len(samples)}")
    print(f"Saved complete subset to {SUBSET_PATH.name}")


if __name__ == "__main__":
    main()
