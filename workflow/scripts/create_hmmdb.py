import os
import argparse
from ftplib import FTP
import subprocess


def download_rfam_files(rfam_ids, outpath):
    """Download Rfam files for the given IDs if they do not already exist."""
    ftp_server = "ftp.ebi.ac.uk"
    ftp_path = "/pub/databases/Rfam/CURRENT/fasta_files/"

    # Connect to the FTP server
    ftp = FTP(ftp_server)
    ftp.login()

    # Download each file for the IDs
    for id in rfam_ids:
        file_name = f"{id}.fa.gz"
        local_file_path = os.path.join(outpath, file_name)

        # Check if the file already exists
        if os.path.exists(local_file_path):
            print(f"File already exists: {file_name}, skipping download.")
            continue

        try:
            with open(local_file_path, "wb") as local_file:
                ftp.retrbinary(f"RETR {ftp_path}{file_name}", local_file.write)
            print(f"Downloaded: {file_name}")
        except Exception as e:
            print(f"Failed to download {file_name}: {e}")

    # Close the FTP connection
    ftp.quit()
    print("Closed FTP connection.")


def run_cd_hit(id, outpath):
    """Run CD-HIT on the downloaded file."""
    input_file = os.path.join(outpath, f"{id}.fa.gz")
    output_file = os.path.join(outpath, f"{id}_97.fasta")

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"Output file already exists: {output_file}, skipping CD-HIT processing.")
        return

    docker_command = [
        "docker",
        "run",
        "-u",
        str(os.getuid()),
        "--rm",
        "-i",
        "-v",
        f"{os.getcwd()}:/data",
        "-w",
        "/data",
        "itvdsbioinfo/cdhit:4.8.1",
        "cd-hit",
        "-c",
        "0.97",
        "-n",
        "5",
        "-T",
        "24",
        "-M",
        "200000",
        "-i",
        input_file,
        "-o",
        output_file,
    ]

    # Execute the Docker command
    try:
        subprocess.run(docker_command, check=True)
        print(f"Processed: {input_file} with CD-HIT, output saved as {output_file}")
    except Exception as e:
        print(f"Failed to process {input_file} with CD-HIT: {e}")


def run_clustalo(id, outpath):
    """Run Clustal Omega on the CD-HIT output file."""
    input_file = os.path.join(outpath, f"{id}_97.fasta")
    output_file = os.path.join(outpath, f"{id}_97.sto")

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(
            f"Output file already exists: {output_file}, skipping Clustal Omega processing."
        )
        return

    docker_command = [
        "docker",
        "run",
        "-u",
        str(os.getuid()),
        "--rm",
        "-i",
        "-v",
        f"{os.getcwd()}:/data",
        "-w",
        "/data",
        "biocontainers/clustalo:v1.2.4-2-deb_cv1",
        "clustalo",
        "-t",
        "DNA",
        "-i",
        input_file,
        "-o",
        output_file,
        "--outfmt=st",
        "--force",
        "--residuenumber",
        "--threads=24",
    ]

    # Execute the Docker command
    try:
        subprocess.run(docker_command, check=True)
        print(
            f"Processed: {input_file} with Clustal Omega, output saved as {output_file}"
        )
    except Exception as e:
        print(f"Failed to process {input_file} with Clustal Omega: {e}")


def run_hmmbuild(id, outpath):
    """Run HMMER's hmmbuild on the Clustal Omega output file."""
    input_file = os.path.join(outpath, f"{id}_97.sto")
    output_file = os.path.join(outpath, f"{id}_97.hmm")

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"Output file already exists: {output_file}, skipping HMMER processing.")
        return

    docker_command = [
        "docker",
        "run",
        "-u",
        str(os.getuid()),
        "--rm",
        "-i",
        "-v",
        f"{os.getcwd()}:/data",
        "-w",
        "/data",
        "staphb/hmmer:3.4",
        "hmmbuild",
        output_file,
        input_file,
    ]

    # Execute the Docker command
    try:
        subprocess.run(docker_command, check=True)
        print(f"Processed: {input_file} with HMMER, output saved as {output_file}")
    except Exception as e:
        print(f"Failed to process {input_file} with HMMER: {e}")


def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Download Rfam files for given IDs.")
    parser.add_argument(
        "--rfam_ids",
        type=str,
        help="Path to the input file containing IDs, one per line.",
    )
    parser.add_argument(
        "--outpath", type=str, help="Path to the folder where files will be saved."
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Initialize an empty list to store the IDs
    id_list = []

    # Read Rfam IDs from file
    try:
        with open(args.rfam_ids, "r") as file:
            for line in file:
                id_list.append(line.strip())
    except FileNotFoundError:
        print(f"Error: The file '{args.rfam_ids}' was not found.")
        exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        exit(1)

    # Create output folder if it doesn't exist
    os.makedirs(args.outpath, exist_ok=True)

    # Download the Rfam files
    download_rfam_files(id_list, args.outpath)

    # Process each file with CD-HIT
    for id in id_list:
        run_cd_hit(id, args.outpath)
        run_clustalo(id, args.outpath)
        run_hmmbuild(id, args.outpath)


if __name__ == "__main__":
    main()
