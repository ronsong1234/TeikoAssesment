# Teiko Immune Cell Analysis

This project loads the supplied clinical-trial cell counts into SQLite,
calculates per-sample cell frequencies, compares miraclib responders with
non-responders, and examines the requested baseline melanoma subset.

## Dashboard

[Open the interactive dashboard](https://ronsong1234.github.io/TeikoAssesment/)

The dashboard includes:

- a sample-level frequency explorer;
- responder versus non-responder boxplots and statistical results; and
- baseline sample counts with a searchable sample table.

## Run in GitHub Codespaces

Open the repository in Codespaces, then run:

```bash
make setup
make pipeline
make dashboard
```

`make dashboard` starts the dashboard server on port 8000. Codespaces will
offer to open the forwarded port in a browser.

The project uses only the Python standard library, so `make setup` does not
download any packages.

## Pipeline

The full pipeline runs these programs in order:

1. `load_data.py` recreates `cell-count.db` and loads `cell-count.csv`.
2. `analysis.py` calculates the relative frequency of every cell population in
   every sample.
3. `statistical_analysis.py` compares melanoma patients receiving miraclib who
   responded with those who did not.
4. `subset_analysis.py` selects baseline melanoma PBMC samples from
   miraclib-treated patients.
5. `generate_dashboard_data.py` prepares the compact JSON file consumed by the
   dashboard.

Every script resolves files relative to its own location, so the pipeline does
not depend on the caller's current directory.

## Output files

| File | Description |
| --- | --- |
| `cell-count.db` | Populated SQLite database |
| `cell_population_frequencies.csv` | Complete Part 2 frequency table |
| `statistical_results.csv` | Part 3 test results and effect sizes |
| `responder_boxplots.svg` | Reproducible responder comparison plot |
| `baseline_melanoma_pbmc_miraclib.csv` | Complete Part 4 baseline subset |
| `dashboard/dashboard-data.json` | Data used by the dashboard |

## Database schema

The database uses five main tables:

```text
projects
  project_id (primary key)

subjects
  subject_id (primary key)
  project_id (foreign key)
  condition, age, sex, treatment, response

samples
  sample_id (primary key)
  subject_id (foreign key)
  sample_type, time_from_treatment_start

cell_types
  cell_type_id (primary key)
  name, display_name

cell_counts
  sample_id (foreign key)
  cell_type_id (foreign key)
  cell_count
  primary key (sample_id, cell_type_id)
```

Project, subject, and sample information are stored once instead of being
repeated for every cell population. Cell counts use a long table, so another
population can be added by inserting a cell type and its counts rather than
altering the database schema.

Foreign keys enforce the project-to-subject-to-sample relationships. Check
constraints reject invalid ages, negative counts, and unsupported response
values. Indexes cover common filters such as project, condition, treatment,
response, sex, subject, timepoint, and cell type.

Two database views simplify analysis:

- `sample_cell_counts` reconstructs the wide shape of the input CSV.
- `cell_population_frequencies` calculates total counts and percentages with a
  window function.

This design works comfortably for hundreds of projects and thousands or
millions of count records because filters and joins operate on indexed keys,
and subject metadata is not duplicated. Loads are performed in one transaction
and the finished database replaces the previous file atomically. For a
multi-user service with heavy concurrent writes, the same logical schema could
be moved to PostgreSQL. At substantially larger analytical scale, the count
table could also be exported to columnar storage or a warehouse while SQLite or
PostgreSQL remains the source of truth.

## Statistical approach

Part 3 includes only PBMC samples from melanoma patients receiving miraclib.
Each patient's relative frequencies are averaged across days 0, 7, and 14
before testing. This avoids treating three repeated samples from one patient as
three independent patients.

Responder and non-responder distributions are compared with a two-sided
Mann-Whitney U test. Benjamini-Hochberg correction controls the false discovery
rate across the five populations, and rank-biserial correlation is reported as
an effect size.

CD4 T cells have an unadjusted `p` value of 0.012, but the adjusted `q` value is
0.062. The difference is therefore nominally significant but does not remain
significant after multiple-testing correction. This is evidence of an
association worth following up, not proof of predictive performance.

## Code structure

The project keeps database loading, descriptive analysis, statistical analysis,
subset analysis, and dashboard-data preparation in separate scripts. This
makes each step easy to run and test on its own while the Makefile still
provides a single reproducible pipeline.

The dashboard is a single `dashboard/index.html` file. It reads generated JSON
rather than connecting directly to SQLite, which keeps the site read-only and
avoids a separate web framework or build step.
