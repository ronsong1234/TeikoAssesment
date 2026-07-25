"""Compare cell frequencies for miraclib responders and non-responders.

The analysis is restricted to PBMC samples from melanoma patients receiving
miraclib. Each patient's percentages are averaged across their samples before
testing so repeated measurements are not treated as independent patients.
"""

import csv
import math
import sqlite3
from pathlib import Path
from statistics import median


ROOT_DIR = Path(__file__).resolve().parent
DATABASE_PATH = ROOT_DIR / "cell-count.db"
RESULTS_PATH = ROOT_DIR / "statistical_results.csv"
LONGITUDINAL_RESULTS_PATH = ROOT_DIR / "longitudinal_results.csv"
BOXPLOT_PATH = ROOT_DIR / "responder_boxplots.svg"

ALPHA = 0.05
POPULATION_ORDER = (
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
)

SUBJECT_FREQUENCY_QUERY = """
SELECT
    sub.subject_id,
    sub.response,
    frequency.population,
    AVG(frequency.percentage) AS average_percentage
FROM cell_population_frequencies AS frequency
JOIN samples AS sample ON sample.sample_id = frequency.sample
JOIN subjects AS sub ON sub.subject_id = sample.subject_id
WHERE sub.condition = 'melanoma'
  AND sub.treatment = 'miraclib'
  AND sub.response IN ('yes', 'no')
  AND sample.sample_type = 'PBMC'
GROUP BY sub.subject_id, sub.response, frequency.population
ORDER BY frequency.population, sub.response, sub.subject_id
"""

LONGITUDINAL_FREQUENCY_QUERY = """
SELECT
    sub.subject_id,
    sub.response,
    sample.time_from_treatment_start,
    frequency.population,
    frequency.percentage
FROM cell_population_frequencies AS frequency
JOIN samples AS sample ON sample.sample_id = frequency.sample
JOIN subjects AS sub ON sub.subject_id = sample.subject_id
WHERE sub.condition = 'melanoma'
  AND sub.treatment = 'miraclib'
  AND sub.response IN ('yes', 'no')
  AND sample.sample_type = 'PBMC'
  AND sample.time_from_treatment_start IN (0, 7, 14)
ORDER BY
    frequency.population,
    sample.time_from_treatment_start,
    sub.response,
    sub.subject_id
"""


def load_subject_frequencies():
    if not DATABASE_PATH.exists():
        raise FileNotFoundError(
            "cell-count.db was not found. Run 'python load_data.py' first."
        )

    connection = sqlite3.connect(DATABASE_PATH)
    try:
        rows = connection.execute(SUBJECT_FREQUENCY_QUERY).fetchall()
    finally:
        connection.close()

    grouped = {
        population: {"yes": [], "no": []}
        for population in POPULATION_ORDER
    }
    for _subject, response, population, percentage in rows:
        grouped[population][response].append(percentage)
    return grouped


def load_longitudinal_frequencies():
    """Return per-timepoint values and within-subject day-14 changes."""
    if not DATABASE_PATH.exists():
        raise FileNotFoundError(
            "cell-count.db was not found. Run 'python load_data.py' first."
        )

    connection = sqlite3.connect(DATABASE_PATH)
    try:
        rows = connection.execute(LONGITUDINAL_FREQUENCY_QUERY).fetchall()
    finally:
        connection.close()

    timepoints = {
        population: {
            timepoint: {"yes": [], "no": []}
            for timepoint in (0, 7, 14)
        }
        for population in POPULATION_ORDER
    }
    subject_values = {}

    for subject, response, timepoint, population, percentage in rows:
        timepoints[population][timepoint][response].append(percentage)
        subject_values.setdefault(
            (subject, response, population), {}
        )[timepoint] = percentage

    changes = {
        population: {"yes": [], "no": []}
        for population in POPULATION_ORDER
    }
    for (_subject, response, population), values in subject_values.items():
        if 0 in values and 14 in values:
            changes[population][response].append(values[14] - values[0])

    return timepoints, changes


def percentile(values, proportion):
    """Calculate a percentile using linear interpolation."""
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot calculate a percentile for an empty group")
    position = (len(ordered) - 1) * proportion
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + (ordered[upper] - ordered[lower]) * fraction


def mann_whitney_u(first_group, second_group):
    """Return U, a two-sided tie-corrected p-value, and rank-biserial effect."""
    n_first = len(first_group)
    n_second = len(second_group)
    if n_first == 0 or n_second == 0:
        raise ValueError("Both comparison groups must contain observations")

    combined = [(value, 0) for value in first_group]
    combined.extend((value, 1) for value in second_group)
    combined.sort(key=lambda item: item[0])

    rank_sum_first = 0.0
    tie_correction = 0
    position = 0
    while position < len(combined):
        tie_end = position + 1
        while (
            tie_end < len(combined)
            and combined[tie_end][0] == combined[position][0]
        ):
            tie_end += 1

        average_rank = ((position + 1) + tie_end) / 2
        for index in range(position, tie_end):
            if combined[index][1] == 0:
                rank_sum_first += average_rank

        tie_size = tie_end - position
        tie_correction += tie_size**3 - tie_size
        position = tie_end

    u_statistic = rank_sum_first - n_first * (n_first + 1) / 2
    total_size = n_first + n_second
    mean_u = n_first * n_second / 2
    variance_u = (
        n_first
        * n_second
        / 12
        * (
            total_size
            + 1
            - tie_correction / (total_size * (total_size - 1))
        )
    )

    if variance_u == 0:
        p_value = 1.0
    else:
        corrected_distance = max(abs(u_statistic - mean_u) - 0.5, 0)
        z_score = corrected_distance / math.sqrt(variance_u)
        p_value = math.erfc(z_score / math.sqrt(2))

    rank_biserial = 2 * u_statistic / (n_first * n_second) - 1
    return u_statistic, p_value, rank_biserial


def adjust_p_values(p_values):
    """Apply the Benjamini-Hochberg false-discovery-rate correction."""
    count = len(p_values)
    ordered_indexes = sorted(range(count), key=lambda index: p_values[index])
    adjusted = [0.0] * count
    running_minimum = 1.0

    for rank in range(count, 0, -1):
        index = ordered_indexes[rank - 1]
        candidate = p_values[index] * count / rank
        running_minimum = min(running_minimum, candidate)
        adjusted[index] = min(running_minimum, 1.0)
    return adjusted


def run_tests(grouped):
    results = []
    for population in POPULATION_ORDER:
        responders = grouped[population]["yes"]
        nonresponders = grouped[population]["no"]
        u_statistic, p_value, effect_size = mann_whitney_u(
            responders, nonresponders
        )
        results.append(
            {
                "population": population,
                "responder_subjects": len(responders),
                "nonresponder_subjects": len(nonresponders),
                "responder_median_pct": median(responders),
                "nonresponder_median_pct": median(nonresponders),
                "median_difference_pct_points": (
                    median(responders) - median(nonresponders)
                ),
                "mann_whitney_u": u_statistic,
                "p_value": p_value,
                "rank_biserial_effect": effect_size,
            }
        )

    adjusted = adjust_p_values([result["p_value"] for result in results])
    for result, q_value in zip(results, adjusted):
        result["significant_p_0_05"] = result["p_value"] < ALPHA
        result["fdr_q_value"] = q_value
        result["significant_fdr_0_05"] = q_value < ALPHA
    return results


def run_timepoint_tests(grouped):
    """Compare responders and non-responders separately at each timepoint."""
    results = []
    for population in POPULATION_ORDER:
        for timepoint in (0, 7, 14):
            responders = grouped[population][timepoint]["yes"]
            nonresponders = grouped[population][timepoint]["no"]
            u_statistic, p_value, effect_size = mann_whitney_u(
                responders, nonresponders
            )
            results.append(
                {
                    "population": population,
                    "timepoint": timepoint,
                    "responder_subjects": len(responders),
                    "nonresponder_subjects": len(nonresponders),
                    "responder_median_pct": median(responders),
                    "nonresponder_median_pct": median(nonresponders),
                    "median_difference_pct_points": (
                        median(responders) - median(nonresponders)
                    ),
                    "mann_whitney_u": u_statistic,
                    "p_value": p_value,
                    "rank_biserial_effect": effect_size,
                }
            )
    return results


def write_results(results):
    fieldnames = (
        "population",
        "responder_subjects",
        "nonresponder_subjects",
        "responder_median_pct",
        "nonresponder_median_pct",
        "median_difference_pct_points",
        "mann_whitney_u",
        "p_value",
        "significant_p_0_05",
        "fdr_q_value",
        "rank_biserial_effect",
        "significant_fdr_0_05",
    )
    with RESULTS_PATH.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)


def write_longitudinal_results(timepoint_results, change_results):
    fieldnames = (
        "analysis",
        "timepoint",
        "population",
        "responder_subjects",
        "nonresponder_subjects",
        "responder_median_pct",
        "nonresponder_median_pct",
        "median_difference_pct_points",
        "mann_whitney_u",
        "p_value",
        "fdr_q_value",
        "rank_biserial_effect",
        "significant_fdr_0_05",
    )
    with LONGITUDINAL_RESULTS_PATH.open(
        "w", newline="", encoding="utf-8"
    ) as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()

        for result in timepoint_results:
            writer.writerow(
                {
                    "analysis": "timepoint",
                    "timepoint": result["timepoint"],
                    "population": result["population"],
                    "responder_subjects": result["responder_subjects"],
                    "nonresponder_subjects": result["nonresponder_subjects"],
                    "responder_median_pct": result["responder_median_pct"],
                    "nonresponder_median_pct": result[
                        "nonresponder_median_pct"
                    ],
                    "median_difference_pct_points": result[
                        "median_difference_pct_points"
                    ],
                    "mann_whitney_u": result["mann_whitney_u"],
                    "p_value": result["p_value"],
                    "rank_biserial_effect": result[
                        "rank_biserial_effect"
                    ],
                }
            )

        for result in change_results:
            writer.writerow(
                {
                    "analysis": "day14_minus_day0",
                    "population": result["population"],
                    "responder_subjects": result["responder_subjects"],
                    "nonresponder_subjects": result[
                        "nonresponder_subjects"
                    ],
                    "responder_median_pct": result[
                        "responder_median_pct"
                    ],
                    "nonresponder_median_pct": result[
                        "nonresponder_median_pct"
                    ],
                    "median_difference_pct_points": result[
                        "median_difference_pct_points"
                    ],
                    "mann_whitney_u": result["mann_whitney_u"],
                    "p_value": result["p_value"],
                    "fdr_q_value": result["fdr_q_value"],
                    "rank_biserial_effect": result[
                        "rank_biserial_effect"
                    ],
                    "significant_fdr_0_05": result[
                        "significant_fdr_0_05"
                    ],
                }
            )


def box_statistics(values):
    q1 = percentile(values, 0.25)
    med = percentile(values, 0.50)
    q3 = percentile(values, 0.75)
    iqr = q3 - q1
    lower_limit = q1 - 1.5 * iqr
    upper_limit = q3 + 1.5 * iqr
    included = [
        value for value in values if lower_limit <= value <= upper_limit
    ]
    outliers = [
        value for value in values if value < lower_limit or value > upper_limit
    ]
    return {
        "q1": q1,
        "median": med,
        "q3": q3,
        "lower_whisker": min(included),
        "upper_whisker": max(included),
        "outliers": outliers,
    }


def write_boxplot(grouped, results):
    """Create a dependency-free SVG boxplot for the five populations."""
    width = 1200
    height = 720
    left = 85
    right = 35
    top = 75
    bottom = 110
    plot_width = width - left - right
    plot_height = height - top - bottom

    all_values = [
        value
        for population in POPULATION_ORDER
        for response in ("yes", "no")
        for value in grouped[population][response]
    ]
    y_max = max(5, math.ceil(max(all_values) * 1.08 / 5) * 5)

    def y_position(value):
        return top + plot_height - value / y_max * plot_height

    result_by_population = {
        result["population"]: result for result in results
    }
    group_width = plot_width / len(POPULATION_ORDER)
    colors = {"yes": "#2878b5", "no": "#e07a2d"}
    labels = {"yes": "Responders", "no": "Non-responders"}
    offsets = {"yes": -24, "no": 24}

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}" '
        'role="img" aria-labelledby="title description">',
        '<title id="title">Cell population frequencies by response</title>',
        '<desc id="description">Boxplots compare subject-average relative '
        "frequencies for miraclib responders and non-responders among melanoma "
        "patients with PBMC samples.</desc>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<g font-family="Arial, sans-serif" fill="#222222">',
        '<text x="85" y="32" font-size="22" font-weight="bold">'
        "PBMC cell frequencies in melanoma patients receiving miraclib</text>",
        '<text x="85" y="55" font-size="13" fill="#555555">'
        "Each observation is one patient’s mean across days 0, 7, and 14; "
        "q-values use Benjamini–Hochberg correction.</text>",
    ]

    tick_step = 5 if y_max <= 40 else 10
    tick = 0
    while tick <= y_max:
        y = y_position(tick)
        svg.extend(
            [
                f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" '
                f'y2="{y:.2f}" stroke="#dddddd" stroke-width="1"/>',
                f'<text x="{left - 12}" y="{y + 4:.2f}" text-anchor="end" '
                f'font-size="12">{tick}</text>',
            ]
        )
        tick += tick_step

    svg.extend(
        [
            f'<line x1="{left}" y1="{top}" x2="{left}" '
            f'y2="{top + plot_height}" stroke="#333333"/>',
            f'<line x1="{left}" y1="{top + plot_height}" '
            f'x2="{width - right}" y2="{top + plot_height}" '
            'stroke="#333333"/>',
            f'<text x="22" y="{top + plot_height / 2}" '
            f'transform="rotate(-90 22 {top + plot_height / 2:.2f})" '
            'text-anchor="middle" font-size="14">'
            "Relative frequency (%)</text>",
        ]
    )

    for population_index, population in enumerate(POPULATION_ORDER):
        center = left + group_width * (population_index + 0.5)
        result = result_by_population[population]
        significance = (
            "FDR significant"
            if result["significant_fdr_0_05"]
            else "FDR ns"
        )
        p_value = result["p_value"]
        q_value = result["fdr_q_value"]
        p_label = "p&lt;0.001" if p_value < 0.001 else f"p={p_value:.3f}"
        q_label = "q&lt;0.001" if q_value < 0.001 else f"q={q_value:.3f}"
        svg.append(
            f'<text x="{center:.2f}" y="{top + 13}" text-anchor="middle" '
            f'font-size="11">{p_label}, {q_label} ({significance})</text>'
        )

        for response in ("yes", "no"):
            values = grouped[population][response]
            stats = box_statistics(values)
            x = center + offsets[response]
            box_width = 34
            q1_y = y_position(stats["q1"])
            median_y = y_position(stats["median"])
            q3_y = y_position(stats["q3"])
            lower_y = y_position(stats["lower_whisker"])
            upper_y = y_position(stats["upper_whisker"])
            color = colors[response]

            svg.extend(
                [
                    f'<line x1="{x:.2f}" y1="{upper_y:.2f}" x2="{x:.2f}" '
                    f'y2="{q3_y:.2f}" stroke="{color}" stroke-width="1.5"/>',
                    f'<line x1="{x:.2f}" y1="{q1_y:.2f}" x2="{x:.2f}" '
                    f'y2="{lower_y:.2f}" stroke="{color}" stroke-width="1.5"/>',
                    f'<line x1="{x - 10:.2f}" y1="{upper_y:.2f}" '
                    f'x2="{x + 10:.2f}" y2="{upper_y:.2f}" '
                    f'stroke="{color}" stroke-width="1.5"/>',
                    f'<line x1="{x - 10:.2f}" y1="{lower_y:.2f}" '
                    f'x2="{x + 10:.2f}" y2="{lower_y:.2f}" '
                    f'stroke="{color}" stroke-width="1.5"/>',
                    f'<rect x="{x - box_width / 2:.2f}" y="{q3_y:.2f}" '
                    f'width="{box_width}" height="{q1_y - q3_y:.2f}" '
                    f'fill="{color}" fill-opacity="0.28" stroke="{color}" '
                    'stroke-width="1.5"/>',
                    f'<line x1="{x - box_width / 2:.2f}" y1="{median_y:.2f}" '
                    f'x2="{x + box_width / 2:.2f}" y2="{median_y:.2f}" '
                    f'stroke="{color}" stroke-width="2.5"/>',
                ]
            )

            for outlier_index, value in enumerate(stats["outliers"]):
                jitter = ((outlier_index % 5) - 2) * 2
                svg.append(
                    f'<circle cx="{x + jitter:.2f}" '
                    f'cy="{y_position(value):.2f}" r="1.8" '
                    f'fill="{color}" fill-opacity="0.55"/>'
                )

        display_name = population.replace("_", " ")
        svg.append(
            f'<text x="{center:.2f}" y="{top + plot_height + 28}" '
            f'text-anchor="middle" font-size="13">{display_name}</text>'
        )

    legend_y = height - 34
    legend_start = width / 2 - 135
    for index, response in enumerate(("yes", "no")):
        x = legend_start + index * 145
        svg.extend(
            [
                f'<rect x="{x:.2f}" y="{legend_y - 11}" width="16" '
                f'height="12" fill="{colors[response]}" fill-opacity="0.35" '
                f'stroke="{colors[response]}"/>',
                f'<text x="{x + 23:.2f}" y="{legend_y}" '
                f'font-size="13">{labels[response]}</text>',
            ]
        )

    svg.extend(["</g>", "</svg>"])
    BOXPLOT_PATH.write_text("\n".join(svg), encoding="utf-8")


def display_results(results):
    print(
        f"{'population':<12} {'responders':>10} {'non-resp.':>10} "
        f"{'median yes':>11} {'median no':>10} {'difference':>11} "
        f"{'p-value':>11} {'p<.05':>7} {'FDR q':>11} {'FDR<.05':>8}"
    )
    print("-" * 121)
    for result in results:
        print(
            f"{result['population']:<12} "
            f"{result['responder_subjects']:>10} "
            f"{result['nonresponder_subjects']:>10} "
            f"{result['responder_median_pct']:>11.3f} "
            f"{result['nonresponder_median_pct']:>10.3f} "
            f"{result['median_difference_pct_points']:>11.3f} "
            f"{result['p_value']:>11.3g} "
            f"{str(result['significant_p_0_05']):>7} "
            f"{result['fdr_q_value']:>11.3g} "
            f"{str(result['significant_fdr_0_05']):>8}"
        )


def display_longitudinal_results(timepoint_results, change_results):
    print("\nPer-timepoint responder versus non-responder p-values")
    print(f"{'population':<12} {'day 0':>10} {'day 7':>10} {'day 14':>10}")
    print("-" * 46)
    by_key = {
        (result["population"], result["timepoint"]): result
        for result in timepoint_results
    }
    for population in POPULATION_ORDER:
        print(
            f"{population:<12} "
            f"{by_key[(population, 0)]['p_value']:>10.3g} "
            f"{by_key[(population, 7)]['p_value']:>10.3g} "
            f"{by_key[(population, 14)]['p_value']:>10.3g}"
        )

    print("\nWithin-subject change from day 0 to day 14")
    print(
        f"{'population':<12} {'median yes':>12} {'median no':>12} "
        f"{'p-value':>11} {'FDR q':>11} {'FDR<.05':>8}"
    )
    print("-" * 71)
    for result in change_results:
        print(
            f"{result['population']:<12} "
            f"{result['responder_median_pct']:>12.3f} "
            f"{result['nonresponder_median_pct']:>12.3f} "
            f"{result['p_value']:>11.3g} "
            f"{result['fdr_q_value']:>11.3g} "
            f"{str(result['significant_fdr_0_05']):>8}"
        )


def main():
    grouped = load_subject_frequencies()
    results = run_tests(grouped)
    timepoints, changes = load_longitudinal_frequencies()
    timepoint_results = run_timepoint_tests(timepoints)
    change_results = run_tests(changes)
    write_results(results)
    write_longitudinal_results(timepoint_results, change_results)
    write_boxplot(grouped, results)
    display_results(results)
    display_longitudinal_results(timepoint_results, change_results)
    print(f"\nSaved detailed results to {RESULTS_PATH.name}")
    print(
        f"Saved longitudinal results to {LONGITUDINAL_RESULTS_PATH.name}"
    )
    print(f"Saved boxplots to {BOXPLOT_PATH.name}")


if __name__ == "__main__":
    main()
