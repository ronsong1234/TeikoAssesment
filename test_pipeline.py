"""Small regression tests for the generated database and analyses."""

import csv
import sqlite3
import unittest
from pathlib import Path

from statistical_analysis import (
    load_longitudinal_frequencies,
    mann_whitney_u,
    run_tests,
)


ROOT_DIR = Path(__file__).resolve().parent
DATABASE_PATH = ROOT_DIR / "cell-count.db"
CSV_PATH = ROOT_DIR / "cell-count.csv"


class PipelineTests(unittest.TestCase):
    def test_database_row_counts_match_input(self):
        with CSV_PATH.open(newline="", encoding="utf-8-sig") as input_file:
            csv_rows = sum(1 for _ in csv.DictReader(input_file))

        with sqlite3.connect(DATABASE_PATH) as connection:
            samples = connection.execute(
                "SELECT COUNT(*) FROM samples"
            ).fetchone()[0]
            counts = connection.execute(
                "SELECT COUNT(*) FROM cell_counts"
            ).fetchone()[0]

        self.assertEqual(samples, csv_rows)
        self.assertEqual(counts, csv_rows * 5)

    def test_sample_percentages_sum_to_100(self):
        with sqlite3.connect(DATABASE_PATH) as connection:
            minimum, maximum = connection.execute(
                """
                SELECT MIN(total_percentage), MAX(total_percentage)
                FROM (
                    SELECT sample, SUM(percentage) AS total_percentage
                    FROM cell_population_frequencies
                    GROUP BY sample
                )
                """
            ).fetchone()

        self.assertAlmostEqual(minimum, 100.0, places=10)
        self.assertAlmostEqual(maximum, 100.0, places=10)

    def test_mann_whitney_textbook_ordering(self):
        statistic, p_value, effect = mann_whitney_u(
            [1, 2, 3], [4, 5, 6]
        )
        self.assertEqual(statistic, 0.0)
        self.assertAlmostEqual(p_value, 0.08085559837005228)
        self.assertEqual(effect, -1.0)

    def test_b_cell_change_survives_fdr_correction(self):
        _timepoints, changes = load_longitudinal_frequencies()
        results = run_tests(changes)
        b_cell = next(
            result for result in results if result["population"] == "b_cell"
        )
        self.assertLess(b_cell["fdr_q_value"], 0.05)
        self.assertTrue(b_cell["significant_fdr_0_05"])

    def test_generated_csv_files_use_lf_line_endings(self):
        for filename in (
            "cell_population_frequencies.csv",
            "statistical_results.csv",
            "longitudinal_results.csv",
            "baseline_melanoma_pbmc_miraclib.csv",
        ):
            self.assertNotIn(b"\r\n", (ROOT_DIR / filename).read_bytes())


if __name__ == "__main__":
    unittest.main()
