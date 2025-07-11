import yaml
import shutil
import logging
import sys


def extract_sample_keys(file_path):
    try:
        with open(file_path, "r") as file:
            data = yaml.safe_load(file)

        samples = data.get("samples", {})
        if isinstance(samples, dict):
            return list(samples.keys())
        else:
            return []

    except FileNotFoundError:
        logging.info(
            f"The Snakemake config file '{file_path}' was not found. Hint: Try running the pipeline without the '-rerun' flag first."
        )
        return []

    except yaml.YAMLError as e:
        logging.error(f"Error: Failed to parse YAML file '{file_path}': {e}")
        return []

    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        return []


if __name__ == "__main__":
    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.info(
        "Rerun flag is set. Removing all your samples files in the workflow/results and resources directory."
    )

    yaml_file = "config/snakemake_config.yaml"
    sample_keys = extract_sample_keys(yaml_file)

    for sample in sample_keys:
        logging.info(f"Removing files for sample: {sample}")
        logging.info(f"Removing workflow/results/{sample} files.")

        try:
            shutil.rmtree(f"workflow/results/{sample}", ignore_errors=True)
        except:
            logging.warning(
                f"Failed to remove workflow/results/{sample}. It may not exist."
            )

        logging.info(f"Removing resources/{sample} files.")

        try:
            shutil.rmtree(f"resources/{sample}", ignore_errors=True)
        except:
            logging.warning(f"Failed to remove resources/{sample}. It may not exist.")
