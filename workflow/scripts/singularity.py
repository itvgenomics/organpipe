import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Process some Docker images.")
parser.add_argument(
    "--sifdir", type=str, required=True, help="Path to where .sif images will be pulled"
)

args = parser.parse_args()

docker_images = [
    "itvdsbioinfo/chloe:1.0",
    "nanozoo/mitos:2.0.3--9b425c9",
    "brunomsilva/novoplasty:4.3.1",
    "ghcr.io/marcelauliano/mitohifi:master",
    "staphb/hmmer:3.4",
    "biocontainers/pilon:v1.23dfsg-1-deb_cv1",
    "itvdsbioinfo/hic_mapping:1.0",
    "itvdsbioinfo/ogdraw:1.1.1",
    "pegi3s/blast:2.13.0",
    "staphb/hmmer:3.4",
    "itvdsbioinfo/pimba_adapterremoval:v2.2.3",
    "brunomsilva/cpgavas2:latest"
]

sif_dir = args.sifdir

os.makedirs(os.path.abspath(sif_dir), exist_ok=True)

for image in docker_images:
    print(f"INFO: Fetching {image} to {sif_dir}.")

    # Construct the singularity pull command with output directory
    output_file = os.path.join(
        os.path.abspath(sif_dir), f"{image.split('/')[-1].split(':')[0]}.sif"
    )

    # Check if the SIF file already exists
    if os.path.exists(output_file):
        print(f"INFO: {output_file} already exists. Skipping pull for {image}.")
        continue

    command = f"singularity pull {output_file} docker://{image}"

    # Execute the command
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"INFO: Successfully pulled {image} to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Error pulling {image}: {e}")
