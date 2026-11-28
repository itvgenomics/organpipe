"""Delete the outputs of the samples of a configuration file, so they are processed again (OrganPipe.sh -rerun).

The samples come from the configuration file of the new run, not from config/snakemake_config.yaml,
which still lists the samples of the previous run.
"""

import argparse
import logging
import os
import re
import shutil
import sys

from check_config import VALID_SAMPLE_NAME
from create_snakemake_config import read_config_rows

# Folders that have one subfolder per sample
SAMPLE_DIRS = ["results", "resources", "workflow/reports"]


def sample_dirs(configfile):
    """Existing output folders of the samples of a configuration file."""
    dirs = []
    for row in read_config_rows(configfile):
        sample = row["sample"]
        # A blank or unexpected name would point to a parent folder, such as results/ itself
        if not re.match(VALID_SAMPLE_NAME, sample):
            logging.warning(f"Not deleting the outputs of '{sample}': it is not a valid sample name.")
            continue
        dirs += [os.path.join(parent, sample) for parent in SAMPLE_DIRS if os.path.isdir(os.path.join(parent, sample))]
    return dirs


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--configfile", required=True, help="Configuration file of the new run (.yaml or .csv)")
    parser.add_argument("--dry-run", action="store_true", help="Only list the folders that would be deleted")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )

    dirs = sample_dirs(args.configfile)
    if not dirs:
        logging.info("-rerun: the samples of the configuration file have no earlier outputs to delete.")

    for path in dirs:
        if args.dry_run:
            logging.info(f"-rerun: {path} would be deleted (dry run).")
        else:
            logging.info(f"-rerun: deleting {path}")
            shutil.rmtree(path)


if __name__ == "__main__":
    main()
