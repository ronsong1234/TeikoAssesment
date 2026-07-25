"use client";

import { useEffect, useMemo, useState } from "react";

type BoxValues = {
  lower: number;
  q1: number;
  median: number;
  q3: number;
  upper: number;
};

type TestResult = {
  population: string;
  responders: number;
  nonresponders: number;
  responder_median: number;
  nonresponder_median: number;
  difference: number;
  p_value: number;
  q_value: number;
  effect: number;
  nominal_significant: boolean;
  fdr_significant: boolean;
};

type DashboardData = {
  overview: {
    projects: number;
    subjects: number;
    samples: number;
    populations: number;
  };
  populations: string[];
  samples: [string, number, [number, number][]][];
  response_analysis: {
    boxes: Record<string, { yes: BoxValues; no: BoxValues }>;
    results: TestResult[];
  };
  baseline: {
    total: number;
    projects: Record<string, number>;
    responses: Record<string, number>;
    sexes: Record<string, number>;
    samples: [string, string, string, string, string][];
  };
};

const displayName = (value: string) =>
  value
    .replaceAll("_", " ")
    .replace(/\b\w/g, (letter) => letter.toUpperCase());

const formatPValue = (value: number) =>
  value < 0.001 ? "<0.001" : value.toFixed(3);

function FrequencyExplorer({ data }: { data: DashboardData }) {
  const [sampleIndex, setSampleIndex] = useState(0);
  const [query, setQuery] = useState("sample00000");
  const sample = data.samples[sampleIndex];

  const findSample = () => {
    const index = data.samples.findIndex(
      ([sampleId]) => sampleId.toLowerCase() === query.trim().toLowerCase(),
    );
    if (index >= 0) setSampleIndex(index);
  };

  const move = (direction: number) => {
    const next = Math.min(
      data.samples.length - 1,
      Math.max(0, sampleIndex + direction),
    );
    setSampleIndex(next);
    setQuery(data.samples[next][0]);
  };

  return (
    <section aria-labelledby="frequency-title">
      <div className="section-heading">
        <div>
          <p className="eyebrow">Part 2</p>
          <h2 id="frequency-title">Cell frequency explorer</h2>
          <p>
            Review counts and relative frequencies for any of the 10,500
            samples.
          </p>
        </div>
        <div className="sample-search">
          <label htmlFor="sample-id">Sample ID</label>
          <div>
            <input
              id="sample-id"
              value={query}
              onChange={(event) => setQuery(event.target.value)}
              onKeyDown={(event) => {
                if (event.key === "Enter") findSample();
              }}
            />
            <button type="button" onClick={findSample}>
              Find
            </button>
          </div>
        </div>
      </div>

      <div className="sample-toolbar">
        <button
          type="button"
          onClick={() => move(-1)}
          disabled={sampleIndex === 0}
        >
          Previous
        </button>
        <div>
          <strong>{sample[0]}</strong>
          <span>{sample[1].toLocaleString()} total cells</span>
        </div>
        <button
          type="button"
          onClick={() => move(1)}
          disabled={sampleIndex === data.samples.length - 1}
        >
          Next
        </button>
      </div>

      <div className="frequency-list">
        {data.populations.map((population, index) => {
          const [count, percentage] = sample[2][index];
          return (
            <div className="frequency-row" key={population}>
              <div className="frequency-label">
                <strong>{displayName(population)}</strong>
                <span>{count.toLocaleString()} cells</span>
              </div>
              <div className="bar-track" aria-hidden="true">
                <span style={{ width: `${percentage}%` }} />
              </div>
              <strong className="frequency-value">
                {percentage.toFixed(2)}%
              </strong>
            </div>
          );
        })}
      </div>
    </section>
  );
}

function ResponseBoxplot({ data }: { data: DashboardData }) {
  const y = (value: number) => 300 - (value / 45) * 240;
  const groupWidth = 168;

  return (
    <div className="chart-wrap">
      <svg
        className="boxplot"
        viewBox="0 0 960 360"
        role="img"
        aria-labelledby="boxplot-title boxplot-description"
      >
        <title id="boxplot-title">
          Cell frequency comparison by treatment response
        </title>
        <desc id="boxplot-description">
          Five paired boxplots compare patient-average PBMC frequencies for
          responders and non-responders.
        </desc>
        {[0, 10, 20, 30, 40].map((tick) => (
          <g key={tick}>
            <line x1="72" y1={y(tick)} x2="930" y2={y(tick)} />
            <text x="58" y={y(tick) + 4} textAnchor="end">
              {tick}
            </text>
          </g>
        ))}
        <text className="axis-title" x="18" y="185" textAnchor="middle">
          Frequency (%)
        </text>

        {data.populations.map((population, index) => {
          const center = 120 + index * groupWidth;
          const result = data.response_analysis.results[index];
          return (
            <g key={population}>
              {(["yes", "no"] as const).map((response) => {
                const box = data.response_analysis.boxes[population][response];
                const x = center + (response === "yes" ? -20 : 20);
                const colorClass =
                  response === "yes" ? "responder-box" : "nonresponder-box";
                return (
                  <g className={colorClass} key={response}>
                    <line x1={x} y1={y(box.upper)} x2={x} y2={y(box.q3)} />
                    <line x1={x} y1={y(box.q1)} x2={x} y2={y(box.lower)} />
                    <line
                      x1={x - 8}
                      y1={y(box.upper)}
                      x2={x + 8}
                      y2={y(box.upper)}
                    />
                    <line
                      x1={x - 8}
                      y1={y(box.lower)}
                      x2={x + 8}
                      y2={y(box.lower)}
                    />
                    <rect
                      x={x - 14}
                      y={y(box.q3)}
                      width="28"
                      height={y(box.q1) - y(box.q3)}
                    />
                    <line
                      className="median"
                      x1={x - 14}
                      y1={y(box.median)}
                      x2={x + 14}
                      y2={y(box.median)}
                    />
                  </g>
                );
              })}
              <text className="chart-label" x={center} y="326" textAnchor="middle">
                {displayName(population)}
              </text>
              <text className="chart-pvalue" x={center} y="40" textAnchor="middle">
                p={formatPValue(result.p_value)} · q=
                {formatPValue(result.q_value)}
              </text>
            </g>
          );
        })}
      </svg>
      <div className="legend" aria-label="Chart legend">
        <span><i className="responder-swatch" />Responders</span>
        <span><i className="nonresponder-swatch" />Non-responders</span>
      </div>
    </div>
  );
}

function ResponseAnalysis({ data }: { data: DashboardData }) {
  return (
    <section aria-labelledby="response-title">
      <div className="section-heading">
        <div>
          <p className="eyebrow">Part 3</p>
          <h2 id="response-title">Response analysis</h2>
          <p>
            Melanoma patients receiving miraclib, restricted to PBMC samples.
            Each observation is a patient mean across days 0, 7, and 14.
          </p>
        </div>
        <div className="finding">
          <span>Primary finding</span>
          <strong>CD4 T cells: p=0.012, q=0.062</strong>
          <small>Nominally significant; not significant after FDR correction</small>
        </div>
      </div>

      <ResponseBoxplot data={data} />

      <div className="table-scroll">
        <table>
          <thead>
            <tr>
              <th>Population</th>
              <th>Responder median</th>
              <th>Non-responder median</th>
              <th>Difference</th>
              <th>p-value</th>
              <th>FDR q-value</th>
              <th>Effect</th>
            </tr>
          </thead>
          <tbody>
            {data.response_analysis.results.map((result) => (
              <tr key={result.population}>
                <td><strong>{displayName(result.population)}</strong></td>
                <td>{result.responder_median.toFixed(2)}%</td>
                <td>{result.nonresponder_median.toFixed(2)}%</td>
                <td className={result.difference > 0 ? "positive" : "negative"}>
                  {result.difference > 0 ? "+" : ""}
                  {result.difference.toFixed(2)} pp
                </td>
                <td>{formatPValue(result.p_value)}</td>
                <td>{formatPValue(result.q_value)}</td>
                <td>{result.effect.toFixed(3)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  );
}

function BaselineAnalysis({ data }: { data: DashboardData }) {
  const [query, setQuery] = useState("");
  const [project, setProject] = useState("all");
  const rows = useMemo(() => {
    const normalized = query.trim().toLowerCase();
    return data.baseline.samples
      .filter(
        ([rowProject, subject, , , sample]) =>
          (project === "all" || rowProject === project) &&
          (!normalized ||
            subject.toLowerCase().includes(normalized) ||
            sample.toLowerCase().includes(normalized)),
      )
      .slice(0, 100);
  }, [data.baseline.samples, project, query]);

  return (
    <section aria-labelledby="baseline-title">
      <div className="section-heading">
        <div>
          <p className="eyebrow">Part 4</p>
          <h2 id="baseline-title">Baseline subset</h2>
          <p>
            Melanoma PBMC samples at day 0 from patients treated with miraclib.
          </p>
        </div>
      </div>

      <div className="summary-grid">
        <article>
          <span>Qualifying samples</span>
          <strong>{data.baseline.total}</strong>
        </article>
        <article>
          <span>Projects</span>
          <strong>
            {Object.entries(data.baseline.projects)
              .map(([name, count]) => `${name}: ${count}`)
              .join(" · ")}
          </strong>
        </article>
        <article>
          <span>Response</span>
          <strong>
            {data.baseline.responses.yes} yes · {data.baseline.responses.no} no
          </strong>
        </article>
        <article>
          <span>Sex</span>
          <strong>
            {data.baseline.sexes.F} female · {data.baseline.sexes.M} male
          </strong>
        </article>
      </div>

      <div className="table-tools">
        <label>
          Search sample or subject
          <input
            value={query}
            onChange={(event) => setQuery(event.target.value)}
            placeholder="sample00000"
          />
        </label>
        <label>
          Project
          <select
            value={project}
            onChange={(event) => setProject(event.target.value)}
          >
            <option value="all">All projects</option>
            {Object.keys(data.baseline.projects).map((name) => (
              <option value={name} key={name}>{name}</option>
            ))}
          </select>
        </label>
      </div>

      <div className="table-scroll">
        <table>
          <thead>
            <tr>
              <th>Sample</th>
              <th>Subject</th>
              <th>Project</th>
              <th>Response</th>
              <th>Sex</th>
            </tr>
          </thead>
          <tbody>
            {rows.map(([rowProject, subject, response, sex, sample]) => (
              <tr key={sample}>
                <td><strong>{sample}</strong></td>
                <td>{subject}</td>
                <td>{rowProject}</td>
                <td>{response === "yes" ? "Responder" : "Non-responder"}</td>
                <td>{sex === "F" ? "Female" : "Male"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p className="table-note">
        Showing up to 100 matching rows. The complete subset is included as a
        CSV output in the repository.
      </p>
    </section>
  );
}

export default function Home() {
  const [data, setData] = useState<DashboardData | null>(null);
  const [activeTab, setActiveTab] = useState("frequencies");
  const [error, setError] = useState("");

  useEffect(() => {
    fetch("/dashboard-data.json")
      .then((response) => {
        if (!response.ok) throw new Error("Dashboard data could not be loaded.");
        return response.json();
      })
      .then(setData)
      .catch((requestError) => setError(requestError.message));
  }, []);

  return (
    <main>
      <header className="site-header">
        <div className="brand">
          <span className="brand-mark">LB</span>
          <div>
            <strong>Loblaw Bio</strong>
            <span>Clinical analytics</span>
          </div>
        </div>
        <div className="status"><i />Analysis complete</div>
      </header>

      <div className="hero">
        <p className="eyebrow">Miraclib clinical trial</p>
        <h1>Immune cell population analysis</h1>
        <p>
          Explore sample composition, compare treatment response, and review
          the baseline melanoma cohort.
        </p>
        {data && (
          <div className="hero-stats">
            <span><strong>{data.overview.samples.toLocaleString()}</strong> samples</span>
            <span><strong>{data.overview.subjects.toLocaleString()}</strong> subjects</span>
            <span><strong>{data.overview.projects}</strong> projects</span>
            <span><strong>{data.overview.populations}</strong> populations</span>
          </div>
        )}
      </div>

      <nav className="tabs" aria-label="Dashboard sections">
        {[
          ["frequencies", "Sample frequencies"],
          ["response", "Response analysis"],
          ["baseline", "Baseline subset"],
        ].map(([value, label]) => (
          <button
            type="button"
            key={value}
            aria-current={activeTab === value ? "page" : undefined}
            onClick={() => setActiveTab(value)}
          >
            {label}
          </button>
        ))}
      </nav>

      <div className="content">
        {!data && !error && <p className="loading">Loading trial data…</p>}
        {error && <p className="error">{error}</p>}
        {data && activeTab === "frequencies" && <FrequencyExplorer data={data} />}
        {data && activeTab === "response" && <ResponseAnalysis data={data} />}
        {data && activeTab === "baseline" && <BaselineAnalysis data={data} />}
      </div>

      <footer>
        Relative frequencies are calculated from the five measured immune cell
        populations in each sample.
      </footer>
    </main>
  );
}
