"""How the values of the config file are read. Shared by check_config.py and create_snakemake_config.py."""

import math


def as_text(value):
    """Return a config value as text, with every kind of blank value as "".

    A blank field can arrive as None (`field:` in YAML), "" (`field: ""`, as written by the
    config generator), NaN (an empty CSV cell) or only spaces. YAML reads unquoted Yes/No
    as booleans, so they are turned back into "Yes" and "No".
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    if isinstance(value, bool):
        return "Yes" if value else "No"
    text = str(value)
    return text if text.strip() else ""


def read_layout(config):
    """Layout of the short reads: "paired" (_R1/_R2 files) or "single" (one file of single-end or merged reads).

    The field is optional, so the configs written before it existed keep their paired reads.
    """
    return as_text(config.get("read_layout")).strip().lower() or "paired"


# Fields that choose the annotation tools of chloroplasts, used when `annotation` is "yes".
CHLORO_ANNOTATOR_FIELDS = ["run_cpgavas2", "run_chloe"]


def chloro_annotator(config, field):
    """Whether a chloroplast annotation tool runs (run_cpgavas2 or run_chloe): "yes", "no" or another value to report.

    The fields are optional and "yes" when blank, so the configs written before they existed keep both tools.
    """
    return as_text(config.get(field)).strip().lower() or "yes"
