"""Display the relative frequency of each cell population in every sample."""

import csv
import sqlite3
from pathlib import Path


DATABASE_PATH = Path(__file__).resolve().parent / "cell-count.db"
OUTPUT_PATH = Path(__file__).resolve().parent / "cell_population_frequencies.csv"

FREQUENCY_QUERY = """
SELECT sample, total_count, population, count, percentage
FROM cell_population_frequencies
ORDER BY sample, population
"""


def get_population_frequencies():
    """Return all cell-population frequency rows from the database."""
    if not DATABASE_PATH.exists():
        raise FileNotFoundError(
            "cell-count.db was not found. Run 'python load_data.py' first."
        )

    connection = sqlite3.connect(DATABASE_PATH)
    connection.row_factory = sqlite3.Row
    try:
        return connection.execute(FREQUENCY_QUERY).fetchall()
    finally:
        connection.close()


def display_frequency_table(rows):
    """Print frequency rows as a readable table."""
    headings = ("sample", "total_count", "population", "count", "percentage")
    print(
        f"{headings[0]:<12} "
        f"{headings[1]:>11} "
        f"{headings[2]:<12} "
        f"{headings[3]:>8} "
        f"{headings[4]:>10}"
    )
    print("-" * 57)

    for row in rows:
        print(
            f"{row['sample']:<12} "
            f"{row['total_count']:>11} "
            f"{row['population']:<12} "
            f"{row['count']:>8} "
            f"{row['percentage']:>9.2f}%"
        )


def write_frequency_table(rows):
    """Save the complete frequency table as a CSV file."""
    columns = ("sample", "total_count", "population", "count", "percentage")
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.writer(output_file, lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(
                (
                    row["sample"],
                    row["total_count"],
                    row["population"],
                    row["count"],
                    f"{row['percentage']:.6f}",
                )
            )


def main():
    rows = get_population_frequencies()
    write_frequency_table(rows)
    display_frequency_table(rows[:20])
    print(f"\nShowing 20 of {len(rows)} rows.")
    print(f"Saved complete table to {OUTPUT_PATH.name}")


if __name__ == "__main__":
    main()
