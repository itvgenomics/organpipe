import os
import subprocess
import argparse
import logging
import sys
import yaml
from naming import sif_name

FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

logging.basicConfig(
    level=logging.INFO,
    stream=sys.stdout,
    format=FORMAT,
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Images that do not stop the run when they cannot be pulled, because only some
# configurations use them. cpgavas2 is only used for chloroplast annotation.
OPTIONAL_IMAGES = {"cpgavas2"}


def pull_command(image, output_file):
    return ["singularity", "pull", "--disable-cache", output_file, f"docker://{image}"]


parser = argparse.ArgumentParser(
    description="Pull the Singularity images listed in workflow/containers.yaml."
)
parser.add_argument(
    "--sifdir", type=str, required=True, help="Path to where .sif images will be pulled"
)
parser.add_argument(
    "--containers",
    type=str,
    default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "containers.yaml"
    ),
    help="YAML file with the images to pull (default: workflow/containers.yaml)",
)

args = parser.parse_args()

with open(args.containers) as handle:
    containers = yaml.safe_load(handle) or {}

sif_dir = os.path.abspath(args.sifdir)

os.makedirs(sif_dir, exist_ok=True)

failed = []

for name, image in containers.items():
    output_file = os.path.join(sif_dir, sif_name(name, image))

    # Check if the SIF file already exists
    if os.path.exists(output_file):
        logging.info(f"{output_file} already exists. Skipping pull for {image}.")
        continue

    logging.info(f"Fetching {image} to {output_file}.")

    # Pull to a temporary name, so an interrupted pull never leaves a .sif behind
    partial_file = f"{output_file}.partial"
    if os.path.exists(partial_file):
        os.remove(partial_file)

    result = subprocess.run(pull_command(image, partial_file))

    if result.returncode != 0 or not os.path.exists(partial_file):
        if os.path.exists(partial_file):
            os.remove(partial_file)
        failed.append((name, image, output_file))
        continue

    os.replace(partial_file, output_file)
    logging.info(f"Successfully pulled {image} to {output_file}")

if failed:
    print(
        f"\nCould not pull {len(failed)} image(s). If this machine has no internet "
        "access, run these commands on another machine and copy the files to the "
        "paths shown:",
        file=sys.stderr,
    )
    for name, image, output_file in failed:
        command = " ".join(pull_command(image, os.path.basename(output_file)))
        print(f"  {command}   # -> {output_file}", file=sys.stderr)

    if any(name not in OPTIONAL_IMAGES for name, _, _ in failed):
        print("ERROR: Required images are missing.", file=sys.stderr)
        sys.exit(1)

    print(
        "WARNING: Continuing without the images above. Jobs that use them will fail.",
        file=sys.stderr,
    )
