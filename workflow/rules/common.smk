import re
import yaml

# Container images, pulled by workflow/scripts/singularity.py.
with open(os.path.join(workflow.basedir, "containers.yaml")) as handle:
    CONTAINERS = yaml.safe_load(handle)

def sif(name):
    """Path of an image in workflow/containers.yaml.

    The file name must match sif_name() in workflow/scripts/singularity.py.
    """
    image = CONTAINERS[name]
    tag = image.rsplit(":", 1)[1] if ":" in image.split("/")[-1] else "latest"
    tag = re.sub(r"[^A-Za-z0-9._-]", "_", tag)
    return f"{config['sif_dir']}/{name}_{tag}.sif"
