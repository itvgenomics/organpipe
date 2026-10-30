#!/usr/bin/env python3
"""
Build an HTML report from the Snakemake benchmark files written by OrganPipe.

Every file under the benchmarks directory is read, assigned to the rule that
produced it (patterns are parsed from the `benchmark:` directives in
workflow/rules/*.smk) and aggregated per stage, rule, sample and k-mer.
Nothing is hard-coded about sample names, seeds, k-mers or which steps ran:
the report reflects whatever is present on disk.

Usage (from the OrganPipe working directory):
    python workflow/scripts/benchmark_report.py -c config/config.csv
    python workflow/scripts/benchmark_report.py            # uses config/snakemake_config.yaml if present
"""

import argparse
import csv
import datetime
import html
import logging
import os
import re
import statistics
import sys
from collections import defaultdict


# ── Constants ──────────────────────────────────────────────────────────────────

METRICS = [
    "s",
    "max_rss",
    "max_vms",
    "max_uss",
    "max_pss",
    "io_in",
    "io_out",
    "mean_load",
    "cpu_time",
]

# Canonical pipeline order of the rule files; it also fixes each stage's color
# so a stage keeps its hue whichever stages were run.
STAGE_ORDER = [
    "trimming",
    "assembly",
    "error_correction",
    "annotation",
    "genbank",
    "images",
    "nhmmer",
]

STAGE_LABELS = {
    "trimming": "Trimming",
    "assembly": "Assembly",
    "error_correction": "Error correction",
    "annotation": "Annotation",
    "genbank": "GenBank",
    "images": "Images",
    "nhmmer": "nhmmer",
    "check_files": "Checks",
    "parse_results": "Results",
    "other": "Other",
}

# Used when a benchmark file matches no rule pattern: the first directory under
# the sample folder usually names the step.
DIR_TO_STAGE = {
    "trimming": "trimming",
    "novoplasty": "assembly",
    "getorganelle": "assembly",
    "mitohifi": "assembly",
    "pilon": "error_correction",
    "mitos2": "annotation",
    "cpgavas2": "annotation",
    "chloe": "annotation",
    "genbanks": "genbank",
    "images": "images",
    "nhmmer": "nhmmer",
}

ASSEMBLERS = [
    ("getorganelle", "GetOrganelle"),
    ("mitohifi", "MitoHiFi"),
    ("_long", "MitoHiFi"),
    ("novoplasty", "NOVOPlasty"),
    ("create_hash", "NOVOPlasty"),
]

WILDCARD_REGEX = {"kmer": r"\d+"}

# Config fields that are never useful in the report.
HIDDEN_CONFIG_FIELDS = {"sample", "seeds"}


# ── Parsing ────────────────────────────────────────────────────────────────────


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate an HTML report from OrganPipe benchmark files"
    )
    parser.add_argument(
        "-b", "--benchmarks", default="benchmarks", help="Benchmarks directory"
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="Run configuration (.csv or .yaml). Defaults to config/snakemake_config.yaml when present",
    )
    parser.add_argument(
        "-r", "--rules", default="workflow/rules", help="Directory with the .smk rule files"
    )
    parser.add_argument(
        "-o", "--output", default="workflow/reports/benchmark.html", help="Output HTML"
    )
    parser.add_argument(
        "--all-samples",
        action="store_true",
        help="Include samples found in the benchmarks directory but absent from the config",
    )
    parser.add_argument(
        "--top", type=int, default=15, help="Number of longest jobs to list"
    )
    return parser.parse_args()


def pattern_to_regex(pattern, bench_dirname):
    """Convert a Snakemake benchmark path pattern into a regex over paths
    relative to the benchmarks directory. Returns None for malformed patterns."""
    prefix = bench_dirname.rstrip("/") + "/"
    if pattern.startswith(prefix):
        pattern = pattern[len(prefix):]
    elif pattern.startswith("benchmarks/"):
        pattern = pattern[len("benchmarks/"):]

    parts = re.split(r"(\{[^{}]*\})", pattern)
    seen = set()
    regex = []
    literal_len = 0
    for part in parts:
        if part.startswith("{") and part.endswith("}"):
            name = part[1:-1].split(",")[0].strip()
            if not name.isidentifier():
                return None, 0
            if name in seen:
                regex.append(f"(?P={name})")
            else:
                seen.add(name)
                regex.append(f"(?P<{name}>{WILDCARD_REGEX.get(name, r'[^/]+?')})")
        else:
            if "{" in part or "}" in part:
                return None, 0
            regex.append(re.escape(part))
            literal_len += len(part)
    return re.compile("^" + "".join(regex) + "$"), literal_len


def load_rule_patterns(rules_dir, bench_dirname):
    """Read `rule NAME:` blocks and their `benchmark:` paths from every .smk file."""
    patterns = []
    if not os.path.isdir(rules_dir):
        logging.warning(f"Rules directory not found: {rules_dir}")
        return patterns

    for fname in sorted(os.listdir(rules_dir)):
        if not fname.endswith(".smk"):
            continue
        stage = fname[: -len(".smk")]
        with open(os.path.join(rules_dir, fname)) as handle:
            lines = handle.read().splitlines()

        rule = None
        for i, line in enumerate(lines):
            match = re.match(r"^\s*(?:rule|checkpoint)\s+(\w+)\s*:", line)
            if match:
                rule = match.group(1)
                continue
            match = re.match(r"^\s*benchmark\s*:\s*(.*)$", line)
            if not (match and rule):
                continue
            text = match.group(1)
            j = i
            while not text.strip() and j + 1 < len(lines):
                j += 1
                text = lines[j]
            quoted = re.search(r"[\"']([^\"']+)[\"']", text)
            if not quoted:
                continue
            regex, literal_len = pattern_to_regex(quoted.group(1), bench_dirname)
            if regex is None:
                logging.warning(
                    f"Skipping malformed benchmark pattern in {fname} (rule {rule}): {quoted.group(1)}"
                )
                continue
            patterns.append(
                {"rule": rule, "stage": stage, "regex": regex, "specificity": literal_len}
            )

    # Most specific patterns first so e.g. "{seed}_kmer{kmer}_x" wins over "{seed}_x".
    patterns.sort(key=lambda p: -p["specificity"])
    return patterns


def read_benchmark(path):
    """Return the mean of each metric over the rows (repeats) of a benchmark file."""
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "s" not in rows[0]:
        return None

    values = {}
    for metric in METRICS:
        numbers = []
        for row in rows:
            try:
                numbers.append(float(row.get(metric, "")))
            except (TypeError, ValueError):
                pass
        values[metric] = statistics.mean(numbers) if numbers else None
    if values["s"] is None:
        return None
    return values


def fallback_rule_name(stem):
    """Best-effort rule name from a file name: the trailing lowercase tokens,
    after dropping sample/seed/k-mer identifiers (e.g. 19_1-COI_150710_run_pilon -> run_pilon)."""
    tokens = stem.split("_")
    kept = []
    for token in reversed(tokens):
        if re.fullmatch(r"kmer\d+", token):
            continue
        if not re.fullmatch(r"[a-z][A-Za-z0-9]*", token):
            break
        kept.append(token)
    return "_".join(reversed(kept)) or stem


def detect_assembler(text):
    text = text.lower()
    for needle, label in ASSEMBLERS:
        if needle in text:
            return label
    return "Shared"


def expand_run(wildcards):
    """Add the seed and k-mer of a NOVOPlasty run ({seed}_kmer{kmer}) to the wildcards."""
    seed, separator, kmer = wildcards.get("run", "").rpartition("_kmer")
    if separator:
        wildcards.update(seed=seed, kmer=kmer)
    return wildcards


def classify(relpath, patterns, organelles, known_samples):
    """Return (sample, rule, stage, wildcards) for a benchmark path."""
    candidates = []
    for pattern in patterns:
        match = pattern["regex"].match(relpath)
        if match:
            if candidates and pattern["specificity"] < candidates[0][0]["specificity"]:
                break
            candidates.append((pattern, match))

    if candidates:
        pattern, match = candidates[0]
        wildcards = expand_run(match.groupdict())
        sample = wildcards.get("sample") or relpath.split("/")[0]
        # Rules that share a benchmark path (e.g. mito/chloro variants): use the organelle.
        organelle = organelles.get(sample, "")
        if len(candidates) > 1 and organelle:
            for cand, cand_match in candidates:
                if organelle in cand["rule"].lower():
                    pattern, wildcards = cand, expand_run(cand_match.groupdict())
                    break
        return sample, pattern["rule"], pattern["stage"], wildcards

    parts = relpath.split("/")
    stem = re.sub(r"\.(benchmark|txt|tsv)$", "", parts[-1])
    if len(parts) > 1:
        sample = parts[0]
        stage = DIR_TO_STAGE.get(parts[1] if len(parts) > 2 else "", "other")
    else:
        prefixes = [s for s in known_samples if stem.startswith(s)]
        sample = max(prefixes, key=len) if prefixes else stem.split(".")[0]
        stage = "other"
    if stem.startswith(sample):
        stem = stem[len(sample):].lstrip("._")
    rule = fallback_rule_name(stem)
    if stage == "other":
        for dirname, dir_stage in DIR_TO_STAGE.items():
            if dirname in rule.lower() or dirname in relpath.lower():
                stage = dir_stage
                break
    wildcards = {}
    kmer = re.search(r"(?:^|_|kmer)(\d{2,3})(?=_)", stem)
    if kmer:
        wildcards["kmer"] = kmer.group(1)
    return sample, rule, stage, wildcards


def stringify(value):
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def load_simple_yaml(path):
    """Flat `key: value` reader used when PyYAML is not installed."""
    data = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#") or line.startswith(" "):
                continue
            key, _, value = line.partition(":")
            value = value.split(" #")[0].strip().strip("\"'")
            data[key.strip()] = value
    return data


def load_config(path):
    """Return ({sample: {field: value}}, field order). Supports the CSV config,
    the user YAML config and the generated config/snakemake_config.yaml."""
    if path.endswith(".csv"):
        with open(path, newline="") as handle:
            first = handle.readline()
            handle.seek(0)
            delimiter = ";" if ";" in first else ","
            reader = csv.DictReader(handle, delimiter=delimiter)
            fields = [f for f in reader.fieldnames or [] if f]
            samples = {}
            for row in reader:
                name = (row.get("sample") or "").strip()
                if name:
                    samples[name] = {f: stringify(row.get(f)) for f in fields}
        return samples, fields

    try:
        import yaml

        with open(path) as handle:
            data = yaml.safe_load(handle) or {}
    except ImportError:
        logging.warning("PyYAML not installed; reading the YAML config as flat key/value pairs")
        data = load_simple_yaml(path)

    if isinstance(data.get("samples"), dict):
        samples = {}
        fields = []
        for name, params in data["samples"].items():
            params = params or {}
            for field in params:
                if field not in fields:
                    fields.append(field)
            samples[str(name)] = {f: stringify(v) for f, v in params.items()}
            if isinstance(params.get("seeds"), (list, tuple)):
                samples[str(name)]["n_seeds"] = str(len(params["seeds"]))
        if any("n_seeds" in s for s in samples.values()):
            fields.append("n_seeds")
        return samples, fields

    # User YAML: one parameter set applied to every listed sample ("all" = every sample found).
    fields = [k for k in data if k != "sample"]
    shared = {f: stringify(data[f]) for f in fields}
    names = [s.strip() for s in stringify(data.get("sample")).split(",") if s.strip()]
    if not names or "all" in [n.lower() for n in names]:
        return {"*": shared}, fields
    return {name: dict(shared) for name in names}, fields


# ── Aggregation ────────────────────────────────────────────────────────────────


def summarize(jobs):
    walls = [j["s"] for j in jobs]
    rss = [j["max_rss"] for j in jobs if j["max_rss"] is not None]
    cpu = [j["cpu_time"] for j in jobs if j["cpu_time"] is not None]
    peak_job = max(jobs, key=lambda j: j["max_rss"] or 0)
    return {
        "jobs": len(jobs),
        "wall": sum(walls),
        "mean": statistics.mean(walls),
        "median": statistics.median(walls),
        "max": max(walls),
        "cpu": sum(cpu),
        "rss_median": statistics.median(rss) if rss else None,
        "rss_peak": max(rss) if rss else None,
        "peak_job": peak_job,
        "io_in": sum(j["io_in"] or 0 for j in jobs),
        "io_out": sum(j["io_out"] or 0 for j in jobs),
        "cores": (sum(cpu) / sum(walls)) if cpu and sum(walls) else None,
    }


def assembly_runs(jobs):
    """Distinct assemblies of one sample: each seed (x k-mer) combination seen in any
    job, plus assembly jobs without seed/k-mer wildcards (e.g. a single GetOrganelle run)."""
    runs = {
        ("seed", j["wildcards"]["seed"], j["wildcards"].get("kmer", ""))
        for j in jobs
        if j["wildcards"].get("seed")
    }
    runs |= {("rule", j["rule"]) for j in jobs if j["stage"] == "assembly" and not j["wildcards"]}
    return len(runs)


def coverage_gaps(jobs):
    """Rules benchmarked for fewer seed x k-mer combinations than the rest of the
    same sample and assembler, which points to missing or failed jobs."""
    gaps = []
    seeded = [j for j in jobs if j["wildcards"].get("seed")]
    # A run without assemblies has no jobs after assembly: only compare rules of the same side
    sides = lambda j: (j["sample"], j["assembler"], j["stage"] == "assembly")
    for (sample, _, _), group in group_by(seeded, sides).items():
        combos = {
            rule: {(j["wildcards"]["seed"], j["wildcards"].get("kmer", "")) for j in rule_jobs}
            for rule, rule_jobs in group_by(group, lambda j: j["rule"]).items()
        }
        expected = set().union(*combos.values())
        for rule, found in sorted(combos.items()):
            if len(found) < len(expected):
                gaps.append((sample, rule, len(found), len(expected)))
    return sorted(gaps, key=lambda g: (natural_key(g[0]), g[1]))


def group_by(jobs, key):
    groups = defaultdict(list)
    for job in jobs:
        groups[key(job)].append(job)
    return groups


def stage_sort_key(stage):
    if stage in STAGE_ORDER:
        return (0, STAGE_ORDER.index(stage), stage)
    return (1, 0, stage)


def natural_key(text):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(text))]


# ── Formatting ─────────────────────────────────────────────────────────────────


def esc(value):
    return html.escape(str(value), quote=True)


def fmt_duration(seconds):
    if seconds is None:
        return "–"
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.1f} s"
    if seconds < 3600:
        minutes, sec = divmod(int(round(seconds)), 60)
        return f"{minutes}m {sec:02d}s"
    hours, rem = divmod(int(round(seconds)), 3600)
    if hours < 48:
        return f"{hours}h {rem // 60:02d}m"
    return f"{seconds / 86400:.1f} d"


def fmt_hours(seconds):
    hours = seconds / 3600
    if hours < 10:
        return f"{hours:,.2f} h"
    return f"{hours:,.1f} h"


def fmt_mem(mb):
    if mb is None:
        return "–"
    if mb < 1024:
        return f"{mb:,.0f} MB"
    return f"{mb / 1024:,.2f} GB"


def fmt_num(value, digits=2):
    if value is None:
        return "–"
    return f"{value:,.{digits}f}"


def stage_label(stage):
    return STAGE_LABELS.get(stage, stage.replace("_", " ").capitalize())


def stage_var(stage, stages):
    """CSS color token for a stage: fixed slot by canonical order, then the
    remaining slots for unknown stages, then neutral."""
    if stage in STAGE_ORDER:
        return f"var(--c{STAGE_ORDER.index(stage) + 1})"
    extras = [s for s in stages if s not in STAGE_ORDER]
    index = len(STAGE_ORDER) + extras.index(stage) + 1
    return f"var(--c{index})" if index <= 8 else "var(--c-other)"


def td(value, sort=None, cls="num"):
    sort_attr = f' data-sort="{esc(sort)}"' if sort is not None else ""
    cls_attr = f' class="{cls}"' if cls else ""
    return f"<td{cls_attr}{sort_attr}>{value}</td>"


# ── HTML building blocks ───────────────────────────────────────────────────────


def build_kpis(jobs, samples, expected_samples):
    total = summarize(jobs)
    longest = max(jobs, key=lambda j: j["s"])
    peak = total["peak_job"]
    sample_walls = [sum(j["s"] for j in group) for group in group_by(jobs, lambda j: j["sample"]).values()]

    sample_value = f"{len(samples)}"
    sample_note = "with benchmark data"
    if expected_samples:
        sample_value = f"{len(samples)} / {len(expected_samples)}"
        sample_note = "benchmarked / configured"

    tiles = [
        ("Samples", sample_value, sample_note),
        ("Jobs", f"{len(jobs):,}", f"{len({j['rule'] for j in jobs})} distinct rules"),
        ("Cumulative runtime", fmt_hours(total["wall"]), "sum of job wall-clock times"),
        ("CPU time", fmt_hours(total["cpu"]), f"avg. {fmt_num(total['cores'])} cores per running job"),
        (
            "Runtime per sample",
            fmt_duration(statistics.median(sample_walls)),
            f"median · range {fmt_duration(min(sample_walls))} – {fmt_duration(max(sample_walls))}",
        ),
        ("Peak memory", fmt_mem(peak["max_rss"]), f"{peak['rule']} · {peak['sample']}"),
        ("Longest job", fmt_duration(longest["s"]), f"{longest['rule']} · {longest['sample']}"),
    ]
    out = ['<section class="kpis">']
    for label, value, note in tiles:
        out.append(
            f'<div class="kpi"><div class="kpi-label">{esc(label)}</div>'
            f'<div class="kpi-value">{esc(value)}</div>'
            f'<div class="kpi-note" title="{esc(note)}">{esc(note)}</div></div>'
        )
    out.append("</section>")
    return "\n".join(out)


def build_config_section(config, samples, fields):
    """Shared settings as a key/value grid; the varying ones go to the sample table."""
    rows = [config[s] for s in samples if s in config]
    if not rows:
        return "", []

    shared, varying = [], []
    for field in fields:
        if field in HIDDEN_CONFIG_FIELDS:
            continue
        values = {r.get(field, "") for r in rows}
        if values == {""}:
            continue
        if len(values) == 1:
            shared.append((field, values.pop()))
        else:
            varying.append(field)

    if not shared:
        return "", varying

    items = "".join(
        f'<div class="kv"><dt>{esc(k)}</dt><dd title="{esc(v)}">{esc(v)}</dd></div>'
        for k, v in shared
    )
    section = f"""
<section>
  <h2>Run configuration</h2>
  <p class="sub">Settings shared by all {len(rows)} configured samples. Settings that differ between samples are listed in the per-sample table.</p>
  <dl class="kv-grid">{items}</dl>
</section>"""
    return section, varying


def build_stage_section(jobs, stages):
    by_stage = group_by(jobs, lambda j: j["stage"])
    total_wall = sum(j["s"] for j in jobs)
    max_wall = max(sum(j["s"] for j in g) for g in by_stage.values())

    bars = []
    rows = []
    for stage in stages:
        stats = summarize(by_stage[stage])
        share = stats["wall"] / total_wall * 100 if total_wall else 0
        width = stats["wall"] / max_wall * 100 if max_wall else 0
        color = stage_var(stage, stages)
        tip = f"{stage_label(stage)}|Runtime {fmt_hours(stats['wall'])} ({share:.1f}%)|{stats['jobs']:,} jobs · median {fmt_duration(stats['median'])}"
        bars.append(
            f'<div class="hbar-row"><div class="hbar-label">{esc(stage_label(stage))}</div>'
            f'<div class="hbar-track"><div class="hbar" tabindex="0" style="width:{width:.2f}%;background:{color}" data-tip="{esc(tip)}"></div>'
            f'<span class="hbar-value">{esc(fmt_hours(stats["wall"]))} · {share:.1f}%</span></div></div>'
        )
        rows.append(
            "<tr>"
            f'<td><span class="swatch" style="background:{color}"></span>{esc(stage_label(stage))}</td>'
            + td(f"{stats['jobs']:,}", stats["jobs"])
            + td(fmt_hours(stats["wall"]), stats["wall"])
            + td(f"{share:.1f}%", share)
            + td(fmt_duration(stats["median"]), stats["median"])
            + td(fmt_duration(stats["max"]), stats["max"])
            + td(fmt_hours(stats["cpu"]), stats["cpu"])
            + td(fmt_mem(stats["rss_peak"]), stats["rss_peak"] or 0)
            + "</tr>"
        )

    return f"""
<section>
  <h2>Runtime by stage</h2>
  <p class="sub">Cumulative wall-clock time of all jobs in each stage. Jobs run in parallel, so this measures compute spent, not elapsed time.</p>
  <div class="chart hbars">{''.join(bars)}</div>
  <div class="table-wrap"><table class="sortable">
    <thead><tr><th>Stage</th><th class="num">Jobs</th><th class="num">Runtime</th><th class="num">Share</th>
    <th class="num">Median job</th><th class="num">Longest job</th><th class="num">CPU time</th><th class="num">Peak memory</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


def build_rule_section(jobs, stages):
    by_rule = group_by(jobs, lambda j: j["rule"])
    total_wall = sum(j["s"] for j in jobs)
    ordered = sorted(by_rule.items(), key=lambda kv: -sum(j["s"] for j in kv[1]))

    rows = []
    for rule, group in ordered:
        stats = summarize(group)
        stage = group[0]["stage"]
        share = stats["wall"] / total_wall * 100 if total_wall else 0
        color = stage_var(stage, stages)
        rows.append(
            "<tr>"
            f'<td class="mono">{esc(rule)}</td>'
            f'<td class="nowrap"><span class="swatch" style="background:{color}"></span>{esc(stage_label(stage))}</td>'
            f'<td>{esc(group[0]["assembler"])}</td>'
            + td(f"{stats['jobs']:,}", stats["jobs"])
            + td(fmt_hours(stats["wall"]), stats["wall"])
            + td(f"{share:.1f}%", share)
            + td(fmt_duration(stats["median"]), stats["median"])
            + td(fmt_duration(stats["max"]), stats["max"])
            + td(fmt_mem(stats["rss_median"]), stats["rss_median"] or 0)
            + td(fmt_mem(stats["rss_peak"]), stats["rss_peak"] or 0)
            + td(fmt_num(stats["cores"]), stats["cores"] or 0)
            + "</tr>"
        )

    return f"""
<section>
  <h2>Runtime by rule</h2>
  <p class="sub">One row per Snakemake rule, sorted by cumulative runtime. Avg. cores is CPU time divided by wall-clock time. Click a column header to sort.</p>
  <div class="table-wrap"><table class="sortable">
    <thead><tr><th>Rule</th><th>Stage</th><th>Assembler</th><th class="num">Jobs</th><th class="num">Runtime</th>
    <th class="num">Share</th><th class="num">Median job</th><th class="num">Longest job</th>
    <th class="num">Median memory</th><th class="num">Peak memory</th><th class="num">Avg. cores</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


def build_sample_section(jobs, stages, samples, config, varying, missing):
    by_sample = group_by(jobs, lambda j: j["sample"])
    walls = {s: sum(j["s"] for j in by_sample[s]) for s in samples}
    max_wall = max(walls.values()) if walls else 1

    legend = "".join(
        f'<span class="legend-item"><span class="swatch" style="background:{stage_var(s, stages)}"></span>{esc(stage_label(s))}</span>'
        for s in stages
    )

    bars = []
    for sample in samples:
        group = by_sample[sample]
        by_stage = group_by(group, lambda j: j["stage"])
        segments = []
        for stage in stages:
            if stage not in by_stage:
                continue
            wall = sum(j["s"] for j in by_stage[stage])
            pct = wall / max_wall * 100
            tip = f"{sample}|{stage_label(stage)}: {fmt_hours(wall)} ({wall / walls[sample] * 100:.1f}%)|{len(by_stage[stage]):,} jobs"
            segments.append(
                f'<span class="seg" tabindex="0" style="width:{pct:.3f}%;background:{stage_var(stage, stages)}" data-tip="{esc(tip)}"></span>'
            )
        bars.append(
            f'<div class="hbar-row" data-sample="{esc(sample.lower())}"><div class="hbar-label mono">{esc(sample)}</div>'
            f'<div class="hbar-track"><div class="stack">{"".join(segments)}</div>'
            f'<span class="hbar-value">{esc(fmt_hours(walls[sample]))}</span></div></div>'
        )

    runs = {s: assembly_runs(by_sample[s]) for s in samples}
    show_runs = any(runs.values())

    head = ["<th>Sample</th>"] + [f"<th>{esc(f)}</th>" for f in varying]
    if show_runs:
        head.append('<th class="num" title="Distinct assemblies benchmarked (seed × k-mer combinations, or single assembler runs)">Assemblies</th>')
    head += [
        '<th class="num">Jobs</th>',
        '<th class="num">Runtime</th>',
        '<th class="num">CPU time</th>',
        '<th class="num">Peak memory</th>',
        "<th>Slowest rule</th>",
    ]

    rows = []
    for sample in samples:
        group = by_sample[sample]
        stats = summarize(group)
        params = config.get(sample, {})
        cells = [f'<td class="mono">{esc(sample)}</td>']
        cells += [f'<td title="{esc(params.get(f, ""))}">{esc(params.get(f, "")) or "–"}</td>' for f in varying]
        if show_runs:
            cells.append(td(f"{runs[sample]:,}", runs[sample]))
        rule_walls = {r: sum(j["s"] for j in g) for r, g in group_by(group, lambda j: j["rule"]).items()}
        slowest = max(rule_walls, key=rule_walls.get)
        cells += [
            td(f"{stats['jobs']:,}", stats["jobs"]),
            td(fmt_hours(stats["wall"]), stats["wall"]),
            td(fmt_hours(stats["cpu"]), stats["cpu"]),
            td(fmt_mem(stats["rss_peak"]), stats["rss_peak"] or 0),
            f'<td><span class="mono">{esc(slowest)}</span> <span class="muted">{esc(fmt_hours(rule_walls[slowest]))}</span></td>',
        ]
        rows.append(f'<tr data-sample="{esc(sample.lower())}">' + "".join(cells) + "</tr>")

    for sample in missing:
        params = config.get(sample, {})
        cells = [f'<td class="mono">{esc(sample)}</td>']
        cells += [f'<td>{esc(params.get(f, "")) or "–"}</td>' for f in varying]
        span = (1 if show_runs else 0) + 5
        cells.append(f'<td colspan="{span}" class="warn">No benchmark files found</td>')
        rows.append(f'<tr data-sample="{esc(sample.lower())}" class="missing">' + "".join(cells) + "</tr>")

    return f"""
<section>
  <h2>Samples</h2>
  <div class="section-head">
    <p class="sub">Cumulative runtime per sample, split by stage. Hover a segment for details.</p>
    <input type="search" id="sample-filter" placeholder="Filter samples" aria-label="Filter samples">
  </div>
  <div class="legend">{legend}</div>
  <div class="chart hbars stacked">{''.join(bars)}</div>
  <div class="table-wrap"><table class="sortable" id="sample-table">
    <thead><tr>{''.join(head)}</tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


def build_kmer_section(jobs):
    kmer_jobs = [j for j in jobs if j["wildcards"].get("kmer")]
    kmers = sorted({j["wildcards"]["kmer"] for j in kmer_jobs}, key=natural_key)
    if len(kmers) < 2:
        return ""

    assembly_rules = sorted(
        {j["rule"] for j in kmer_jobs if j["stage"] == "assembly"},
        key=lambda r: -sum(j["s"] for j in kmer_jobs if j["rule"] == r),
    )
    by_kmer = group_by(kmer_jobs, lambda j: j["wildcards"]["kmer"])

    head = ['<th class="num">k-mer</th>', '<th class="num">Jobs</th>', '<th class="num">Runtime</th>']
    head += [f'<th class="num">Median <span class="mono">{esc(r)}</span></th>' for r in assembly_rules]
    head += [f'<th class="num">Peak memory <span class="mono">{esc(r)}</span></th>' for r in assembly_rules[:1]]

    rows = []
    for kmer in kmers:
        group = by_kmer[kmer]
        stats = summarize(group)
        cells = [td(esc(kmer), int(kmer) if kmer.isdigit() else kmer), td(f"{stats['jobs']:,}", stats["jobs"]), td(fmt_hours(stats["wall"]), stats["wall"])]
        for rule in assembly_rules:
            walls = [j["s"] for j in group if j["rule"] == rule]
            cells.append(td(fmt_duration(statistics.median(walls)) if walls else "–", statistics.median(walls) if walls else 0))
        for rule in assembly_rules[:1]:
            rss = [j["max_rss"] for j in group if j["rule"] == rule and j["max_rss"] is not None]
            cells.append(td(fmt_mem(max(rss)) if rss else "–", max(rss) if rss else 0))
        rows.append("<tr>" + "".join(cells) + "</tr>")

    return f"""
<section>
  <h2>Runtime by k-mer</h2>
  <p class="sub">All jobs that carry a k-mer, summed per k-mer size across samples and seeds.</p>
  <div class="table-wrap"><table class="sortable">
    <thead><tr>{''.join(head)}</tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


def build_top_jobs_section(jobs, stages, top):
    ordered = sorted(jobs, key=lambda j: -j["s"])[:top]
    if not ordered:
        return ""
    wildcard_names = sorted(
        {k for j in ordered for k in j["wildcards"] if k != "sample"},
        key=lambda k: (k != "kmer", k),
    )
    head = ["<th>Sample</th>", "<th>Rule</th>", "<th>Stage</th>"]
    head += [f"<th>{esc('k-mer' if w == 'kmer' else w)}</th>" for w in wildcard_names]
    head += ['<th class="num">Wall time</th>', '<th class="num">CPU time</th>', '<th class="num">Peak memory</th>']

    rows = []
    for job in ordered:
        cells = [
            f'<td class="mono">{esc(job["sample"])}</td>',
            f'<td class="mono">{esc(job["rule"])}</td>',
            f'<td class="nowrap"><span class="swatch" style="background:{stage_var(job["stage"], stages)}"></span>{esc(stage_label(job["stage"]))}</td>',
        ]
        cells += [f'<td class="mono">{esc(job["wildcards"].get(w, "")) or "–"}</td>' for w in wildcard_names]
        cells += [
            td(fmt_duration(job["s"]), job["s"]),
            td(fmt_duration(job["cpu_time"]), job["cpu_time"] or 0),
            td(fmt_mem(job["max_rss"]), job["max_rss"] or 0),
        ]
        rows.append("<tr>" + "".join(cells) + "</tr>")

    return f"""
<section>
  <h2>Longest jobs</h2>
  <p class="sub">The {len(ordered)} individual jobs with the highest wall-clock time.</p>
  <div class="table-wrap"><table class="sortable">
    <thead><tr>{''.join(head)}</tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


CSS = """
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10); --hover: rgba(11,11,11,0.035);
  --accent: #2a78d6; --warn: #b54708;
  --c1: #2a78d6; --c2: #eb6834; --c3: #1baf7a; --c4: #eda100;
  --c5: #e87ba4; --c6: #008300; --c7: #4a3aa7; --c8: #e34948; --c-other: #898781;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --hover: rgba(255,255,255,0.04);
    --accent: #3987e5; --warn: #f79009;
    --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500;
    --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --hover: rgba(255,255,255,0.04);
  --accent: #3987e5; --warn: #f79009;
  --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500;
  --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; background: var(--page); color: var(--ink);
  font: 14px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
}
main { max-width: 1200px; margin: 0 auto; padding: 32px 24px 64px; }
header.top { display: flex; justify-content: space-between; align-items: flex-end; gap: 24px; flex-wrap: wrap;
  padding-bottom: 20px; border-bottom: 1px solid var(--grid); margin-bottom: 24px; }
.eyebrow { font-size: 12px; letter-spacing: .08em; text-transform: uppercase; color: var(--muted); margin: 0 0 4px; }
h1 { font-size: 26px; font-weight: 600; margin: 0; letter-spacing: -0.01em; }
h2 { font-size: 17px; font-weight: 600; margin: 0 0 4px; }
.meta { color: var(--ink-2); font-size: 13px; text-align: right; }
.meta { min-width: 0; max-width: 100%; }
.meta div { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 460px; }
section { min-width: 0; background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 20px 22px; margin-bottom: 20px; }
.sub { color: var(--ink-2); margin: 0 0 16px; font-size: 13px; }
.muted { color: var(--muted); }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12.5px; }
.nowrap { white-space: nowrap; }

.kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(150px, 100%), 1fr)); gap: 0; padding: 0; overflow: hidden; }
.kpi { padding: 16px 18px; border-right: 1px solid var(--grid); border-bottom: 1px solid var(--grid); margin: 0 -1px -1px 0; }
.kpi-label { font-size: 12px; color: var(--ink-2); }
.kpi-value { font-size: 22px; font-weight: 600; margin: 2px 0; }
.kpi-note { font-size: 12px; color: var(--muted); overflow-wrap: anywhere; }

.kv-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(min(230px, 100%), 1fr)); gap: 0 24px; margin: 0; }
.kv { display: flex; justify-content: space-between; gap: 12px; padding: 6px 0; border-bottom: 1px solid var(--grid); min-width: 0; }
.kv dt { color: var(--ink-2); font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12px; flex-shrink: 0; }
.kv dd { margin: 0; text-align: right; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; min-width: 0; }

.chart { margin: 4px 0 20px; }
.hbar-row { display: grid; grid-template-columns: 140px 1fr; align-items: center; gap: 12px; min-height: 26px; }
.stacked .hbar-row { min-height: 22px; }
.hbar-label { color: var(--ink-2); font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; text-align: right; }
.hbar-track { display: flex; align-items: center; gap: 8px; min-width: 0; border-left: 1px solid var(--axis); padding: 3px 0; }
.hbar { height: 16px; min-width: 2px; border-radius: 0 4px 4px 0; flex-shrink: 0; }
.stack { display: flex; gap: 2px; flex: 0 1 auto; width: calc(100% - 72px); }
.seg { height: 14px; min-width: 1px; display: block; }
.seg:last-child { border-radius: 0 4px 4px 0; }
.hbar-value { font-size: 12px; color: var(--ink-2); white-space: nowrap; font-variant-numeric: tabular-nums; }
.hbars .hbar-track > .hbar { max-width: calc(100% - 120px); }
.hbar:hover, .seg:hover, .hbar:focus, .seg:focus { filter: brightness(1.12); outline: none; }

.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; margin: 0 0 10px; font-size: 12px; color: var(--ink-2); }
.legend-item { display: inline-flex; align-items: center; }
.swatch { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; flex-shrink: 0; }

.section-head { display: flex; justify-content: space-between; align-items: center; gap: 16px; flex-wrap: wrap; margin-bottom: 12px; }
.section-head .sub { margin: 0; }
input[type=search] { font: inherit; font-size: 13px; padding: 6px 10px; border: 1px solid var(--axis); border-radius: 6px;
  background: var(--page); color: var(--ink); min-width: 200px; }
input[type=search]:focus { outline: 2px solid var(--accent); outline-offset: -1px; }

.table-wrap { overflow-x: auto; border: 1px solid var(--grid); border-radius: 8px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { padding: 7px 10px; text-align: left; border-bottom: 1px solid var(--grid); vertical-align: middle; }
td { max-width: 260px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
tbody tr:last-child td { border-bottom: 0; }
th { font-weight: 600; font-size: 12px; color: var(--ink-2); background: var(--page); position: sticky; top: 0; white-space: nowrap; }
th .mono { font-size: 11.5px; font-weight: 500; }
.sortable th { cursor: pointer; user-select: none; }
.sortable th:hover { color: var(--ink); }
th[aria-sort=ascending]::after { content: " ↑"; }
th[aria-sort=descending]::after { content: " ↓"; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
tbody tr:hover { background: var(--hover); }
.warn { color: var(--warn); }
tr.missing td { color: var(--muted); }
.notes { color: var(--ink-2); font-size: 12.5px; margin: 0; padding-left: 18px; }
.notes li { margin: 3px 0; }
footer { color: var(--muted); font-size: 12px; text-align: center; margin-top: 28px; }

#tip { position: fixed; pointer-events: none; z-index: 10; background: var(--surface); color: var(--ink);
  border: 1px solid var(--border); border-radius: 6px; padding: 8px 10px; font-size: 12px; line-height: 1.45;
  box-shadow: 0 4px 16px rgba(0,0,0,.12); opacity: 0; transition: opacity .08s; max-width: 280px; }
#tip b { display: block; font-weight: 600; margin-bottom: 2px; }
#tip span { display: block; color: var(--ink-2); }

@media (max-width: 640px) {
  main { padding: 20px 16px 48px; }
  section { padding: 16px; }
  .meta { text-align: left; }
  .hbar-row { grid-template-columns: 96px 1fr; gap: 8px; }
}
@media print {
  body { background: #fff; }
  section { break-inside: avoid; border-color: #ddd; }
  input[type=search], #tip { display: none; }
}
"""

JS = """
(function () {
  // Tooltip for chart marks
  var tip = document.getElementById('tip');
  function show(el, x, y) {
    var parts = el.getAttribute('data-tip').split('|');
    tip.innerHTML = '';
    parts.forEach(function (p, i) {
      var n = document.createElement(i === 0 ? 'b' : 'span');
      n.textContent = p; tip.appendChild(n);
    });
    tip.style.opacity = 1;
    var r = tip.getBoundingClientRect();
    var left = Math.min(x + 14, window.innerWidth - r.width - 8);
    var top = y + 14 + r.height > window.innerHeight ? y - r.height - 10 : y + 14;
    tip.style.left = Math.max(8, left) + 'px'; tip.style.top = top + 'px';
  }
  document.querySelectorAll('[data-tip]').forEach(function (el) {
    el.addEventListener('mousemove', function (e) { show(el, e.clientX, e.clientY); });
    el.addEventListener('mouseleave', function () { tip.style.opacity = 0; });
    el.addEventListener('focus', function () { var r = el.getBoundingClientRect(); show(el, r.left + r.width / 2, r.bottom); });
    el.addEventListener('blur', function () { tip.style.opacity = 0; });
  });

  // Sortable tables (uses data-sort when present)
  document.querySelectorAll('table.sortable').forEach(function (table) {
    table.querySelectorAll('thead th').forEach(function (th, col) {
      th.addEventListener('click', function () {
        var asc = th.getAttribute('aria-sort') !== 'ascending';
        table.querySelectorAll('thead th').forEach(function (h) { h.removeAttribute('aria-sort'); });
        th.setAttribute('aria-sort', asc ? 'ascending' : 'descending');
        var body = table.tBodies[0];
        var rows = Array.prototype.slice.call(body.rows);
        function key(row) {
          var cell = row.cells[col];
          if (!cell) return '';
          var v = cell.getAttribute('data-sort');
          if (v === null) v = cell.textContent.trim();
          var n = parseFloat(v);
          return isNaN(n) || !/^-?[\\d.]+(e[-+]?\\d+)?$/i.test(v) ? v.toLowerCase() : n;
        }
        rows.sort(function (a, b) {
          var x = key(a), y = key(b);
          if (typeof x === 'number' && typeof y === 'number') return asc ? x - y : y - x;
          return asc ? String(x).localeCompare(String(y), undefined, {numeric: true})
                     : String(y).localeCompare(String(x), undefined, {numeric: true});
        });
        rows.forEach(function (r) { body.appendChild(r); });
      });
    });
  });

  // Sample filter (chart rows and table rows)
  var filter = document.getElementById('sample-filter');
  if (filter) filter.addEventListener('input', function () {
    var q = filter.value.trim().toLowerCase();
    document.querySelectorAll('[data-sample]').forEach(function (el) {
      el.style.display = !q || el.getAttribute('data-sample').indexOf(q) !== -1 ? '' : 'none';
    });
  });
})();
"""


def build_html(context):
    sections = "\n".join(s for s in context["sections"] if s)
    notes = ""
    if context["notes"]:
        items = "".join(f"<li>{n}</li>" for n in context["notes"])
        notes = f'<section><h2>Notes</h2><ul class="notes">{items}</ul></section>'

    meta = [f"Generated {esc(context['generated'])}", f"Benchmarks: <span class=\"mono\">{esc(context['bench_dir'])}</span>"]
    if context["config_path"]:
        meta.append(f"Config: <span class=\"mono\" title=\"{esc(context['config_path'])}\">{esc(context['config_path'])}</span>")

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OrganPipe Benchmark Report</title>
<style>{CSS}</style>
</head>
<body>
<main>
  <header class="top">
    <div>
      <p class="eyebrow">OrganPipe</p>
      <h1>Benchmark report</h1>
    </div>
    <div class="meta">{''.join(f'<div>{m}</div>' for m in meta)}</div>
  </header>
  {sections}
  {notes}
  <footer>Generated by workflow/scripts/benchmark_report.py from Snakemake benchmark files.</footer>
</main>
<div id="tip" role="tooltip"></div>
<script>{JS}</script>
</body>
</html>
"""


# ── Main ───────────────────────────────────────────────────────────────────────


def main():
    args = parse_arguments()
    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    bench_dir = os.path.normpath(args.benchmarks)
    if not os.path.isdir(bench_dir):
        logging.error(f"Benchmarks directory not found: {bench_dir}")
        sys.exit(1)

    config_path = args.config
    if config_path is None and os.path.exists("config/snakemake_config.yaml"):
        config_path = "config/snakemake_config.yaml"

    config, fields = {}, []
    if config_path:
        if not os.path.exists(config_path):
            logging.error(f"Config file not found: {config_path}")
            sys.exit(1)
        config, fields = load_config(config_path)
        logging.info(f"Loaded config {config_path} ({len(config)} sample entries)")

    patterns = load_rule_patterns(args.rules, os.path.basename(bench_dir))
    logging.info(f"Parsed {len(patterns)} benchmark patterns from {args.rules}")

    organelles = {s: p.get("organelle", "").lower() for s, p in config.items()}
    known_samples = [s for s in config if s != "*"]

    jobs, unreadable, unmatched = [], [], 0
    for root, _, files in os.walk(bench_dir):
        for fname in files:
            path = os.path.join(root, fname)
            relpath = os.path.relpath(path, bench_dir).replace(os.sep, "/")
            try:
                values = read_benchmark(path)
            except (OSError, UnicodeDecodeError, csv.Error):
                values = None
            if values is None:
                unreadable.append(relpath)
                continue
            sample, rule, stage, wildcards = classify(relpath, patterns, organelles, known_samples)
            if not any(p["regex"].match(relpath) for p in patterns):
                unmatched += 1
            wildcards.pop("sample", None)
            jobs.append(
                {
                    **values,
                    "sample": sample,
                    "rule": rule,
                    "stage": stage,
                    "wildcards": wildcards,
                    "assembler": detect_assembler(rule + " " + relpath),
                    "path": relpath,
                }
            )

    # "*" means the YAML config applies to every sample ("sample: all").
    if "*" in config:
        shared = config.pop("*")
        for sample in {j["sample"] for j in jobs}:
            config[sample] = dict(shared)

    expected = list(config) if config else []
    notes = []
    if expected and not args.all_samples:
        extra = sorted({j["sample"] for j in jobs} - set(expected), key=natural_key)
        if extra:
            notes.append(
                f"{len(extra)} sample(s) present in the benchmarks directory but not in the config were excluded "
                f"(use <span class=\"mono\">--all-samples</span> to include them): "
                + ", ".join(f'<span class="mono">{esc(s)}</span>' for s in extra[:20])
                + (" …" if len(extra) > 20 else "")
            )
        jobs = [j for j in jobs if j["sample"] in config]

    if not jobs:
        logging.error("No benchmark data found for the selected samples")
        sys.exit(1)

    samples = sorted({j["sample"] for j in jobs}, key=natural_key)
    missing = [s for s in expected if s not in samples]
    stages = sorted({j["stage"] for j in jobs}, key=stage_sort_key)

    if missing:
        notes.append(f"{len(missing)} configured sample(s) have no benchmark files (failed or not yet run).")
    gaps = coverage_gaps(jobs)
    if gaps:
        shown = "; ".join(
            f'<span class="mono">{esc(s)}</span> · <span class="mono">{esc(r)}</span> {found} of {total}'
            for s, r, found, total in gaps[:15]
        )
        notes.append(
            f"{len(gaps)} rule(s) have fewer benchmarked seed × k-mer combinations than other rules of the same "
            f"sample, suggesting failed, interrupted or removed jobs: {shown}" + (" …" if len(gaps) > 15 else "") + "."
        )
    if unmatched:
        notes.append(
            f"{unmatched} benchmark file(s) did not match any rule pattern in "
            f"<span class=\"mono\">{esc(args.rules)}</span>; their rule names were inferred from the file names."
        )
    if unreadable:
        notes.append(
            f"{len(unreadable)} file(s) could not be read as Snakemake benchmarks and were skipped: "
            + ", ".join(f'<span class="mono">{esc(p)}</span>' for p in unreadable[:10])
            + (" …" if len(unreadable) > 10 else "")
        )
    notes.append(
        "Runtimes are summed over jobs; because Snakemake runs jobs in parallel, they reflect compute "
        "spent rather than elapsed time. Memory is the maximum resident set size (max_rss) of a single job."
    )

    config_section, varying = build_config_section(config, samples + missing, fields)

    context = {
        "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "bench_dir": bench_dir,
        "config_path": config_path,
        "notes": notes,
        "sections": [
            build_kpis(jobs, samples, expected),
            config_section,
            build_stage_section(jobs, stages),
            build_rule_section(jobs, stages),
            build_sample_section(jobs, stages, samples, config, varying, missing),
            build_kmer_section(jobs),
            build_top_jobs_section(jobs, stages, args.top),
        ],
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write(build_html(context))

    logging.info(
        f"Report written to {args.output} ({len(jobs)} jobs, {len(samples)} samples, {len(stages)} stages)"
    )


if __name__ == "__main__":
    main()
