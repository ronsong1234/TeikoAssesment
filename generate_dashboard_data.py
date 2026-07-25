"""Build the JSON file used by the interactive dashboard."""

import json
import sqlite3
from pathlib import Path

from statistical_analysis import (
    POPULATION_ORDER,
    box_statistics,
    load_subject_frequencies,
    run_tests,
)
from subset_analysis import query_database


ROOT_DIR = Path(__file__).resolve().parent
DATABASE_PATH = ROOT_DIR / "cell-count.db"
OUTPUT_PATH = ROOT_DIR / "dashboard" / "public" / "dashboard-data.json"


def get_overview(connection):
    return {
        "projects": connection.execute(
            "SELECT COUNT(*) FROM projects"
        ).fetchone()[0],
        "subjects": connection.execute(
            "SELECT COUNT(*) FROM subjects"
        ).fetchone()[0],
        "samples": connection.execute(
            "SELECT COUNT(*) FROM samples"
        ).fetchone()[0],
        "populations": connection.execute(
            "SELECT COUNT(*) FROM cell_types"
        ).fetchone()[0],
    }


def get_sample_frequencies(connection):
    rows = connection.execute(
        """
        SELECT sample, total_count, population, count, percentage
        FROM cell_population_frequencies
        ORDER BY sample, population
        """
    ).fetchall()

    samples = []
    current_sample = None
    current_total = 0
    values = {}

    for sample, total_count, population, count, percentage in rows:
        if current_sample is not None and sample != current_sample:
            samples.append(
                [
                    current_sample,
                    current_total,
                    [
                        [values[name][0], round(values[name][1], 4)]
                        for name in POPULATION_ORDER
                    ],
                ]
            )
            values = {}

        current_sample = sample
        current_total = total_count
        values[population] = (count, percentage)

    if current_sample is not None:
        samples.append(
            [
                current_sample,
                current_total,
                [
                    [values[name][0], round(values[name][1], 4)]
                    for name in POPULATION_ORDER
                ],
            ]
        )
    return samples


def get_response_analysis():
    grouped = load_subject_frequencies()
    results = run_tests(grouped)

    boxes = {}
    for population in POPULATION_ORDER:
        boxes[population] = {}
        for response in ("yes", "no"):
            stats = box_statistics(grouped[population][response])
            boxes[population][response] = {
                "lower": round(stats["lower_whisker"], 4),
                "q1": round(stats["q1"], 4),
                "median": round(stats["median"], 4),
                "q3": round(stats["q3"], 4),
                "upper": round(stats["upper_whisker"], 4),
            }

    compact_results = []
    for result in results:
        compact_results.append(
            {
                "population": result["population"],
                "responders": result["responder_subjects"],
                "nonresponders": result["nonresponder_subjects"],
                "responder_median": round(
                    result["responder_median_pct"], 4
                ),
                "nonresponder_median": round(
                    result["nonresponder_median_pct"], 4
                ),
                "difference": round(
                    result["median_difference_pct_points"], 4
                ),
                "p_value": result["p_value"],
                "q_value": result["fdr_q_value"],
                "effect": round(result["rank_biserial_effect"], 4),
                "nominal_significant": result["significant_p_0_05"],
                "fdr_significant": result["significant_fdr_0_05"],
            }
        )
    return {"boxes": boxes, "results": compact_results}


def get_baseline_subset():
    samples, projects, responses, sexes = query_database()
    return {
        "total": len(samples),
        "projects": {
            row["project"]: row["sample_count"] for row in projects
        },
        "responses": {
            row["response"]: row["subject_count"] for row in responses
        },
        "sexes": {row["sex"]: row["subject_count"] for row in sexes},
        "samples": [
            [
                row["project"],
                row["subject"],
                row["response"],
                row["sex"],
                row["sample"],
            ]
            for row in samples
        ],
    }


def main():
    if not DATABASE_PATH.exists():
        raise FileNotFoundError(
            "cell-count.db was not found. Run 'python load_data.py' first."
        )

    connection = sqlite3.connect(DATABASE_PATH)
    try:
        data = {
            "overview": get_overview(connection),
            "populations": list(POPULATION_ORDER),
            "samples": get_sample_frequencies(connection),
            "response_analysis": get_response_analysis(),
            "baseline": get_baseline_subset(),
        }
    finally:
        connection.close()

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(
        json.dumps(data, separators=(",", ":")), encoding="utf-8"
    )
    print(f"Saved dashboard data to {OUTPUT_PATH.relative_to(ROOT_DIR)}")


if __name__ == "__main__":
    main()
