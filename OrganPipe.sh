#!/bin/bash
set -e

VERSION="1.3.0"

print_banner() {
cat << EOF

### OrganPipe: An automated tool to facilitate the assembly, annotation, and curation of mitochondrial and chloroplast genomes
#Version: $VERSION
#Authors: Moreira-Oliveira, R. R., Silva, B. M., Molina, M., Oliveira-Lima, M., Leão, T. F., Vasconcelos, S., & Nunes, G. L. (2025).
#https://doi.org/10.21203/RS.3.RS-5686696/V1

###    Copyright (C) 2024  Renato Oliveira
###
###    This program is free software: you can redistribute it and/or modify
###    it under the terms of the GNU General Public License as published by
###    the Free Software Foundation, either version 3 of the License, or
###    any later version.
###
###    This program is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###    GNU General Public License for more details.
###
###    You should have received a copy of the GNU General Public License
###    along with this program.  If not, see <http://www.gnu.org/licenses/>.

EOF
}

usage() {
cat << EOF
OrganPipe $VERSION - assembly, annotation and curation of mitochondrial and chloroplast genomes

Usage:
  bash OrganPipe.sh -d <workdir> -c <configfile> -t <cores> [options]
  bash OrganPipe.sh -d <workdir> -c <configfile> -slurm -j <jobs> -partition <name> [options]

Required:
  -d <path>           Working directory (the OrganPipe directory)
  -c <file>           Configuration file (.yaml or .csv)
  -t <int>            Local mode: number of CPU cores Snakemake may use at once
  -j <int>            SLURM mode: maximum number of jobs queued or running at once

SLURM:
  -slurm              Submit each job to SLURM (resources in config/resources.yaml)
  -partition <name>   Default SLURM partition (required with -slurm)

Options:
  -np                 Dry run: show the jobs that would run without running them
  -unlock             Unlock the working directory after an interrupted run
  -rerun              Delete results/, resources/ and workflow/reports/ of the samples
                      from the previous run before starting
  -batch              Split the DAG into batches (useful for many samples)
  -nbatch <int>       Number of batches used with -batch (default: 3)
  -sifdir <path>      Directory for the Singularity images (default: resources/sif_dir)
  -nhmmer_db <path>   HMM database used by HMMER (default: resources/rfam.hmm)
  -mitos_db <path>    MITOS2 reference directory, inside the working directory
                      (default: resources/refseq89m)
  -mem <int>          Local mode: memory (GB) that NOVOPlasty jobs may use at once
  -h, --help          Show this help and exit
  -v, --version       Show the OrganPipe version and exit

Examples:
  bash OrganPipe.sh -d . -c config/config.csv -t 8 -np
  bash OrganPipe.sh -d . -c config/config.csv -slurm -j 50 -partition <your_partition>

Documentation: https://github.com/itvgenomics/organpipe
EOF
}

fail() {
    echo "ERROR: $1" >&2
    echo "Run 'bash OrganPipe.sh -h' for usage." >&2
    exit 1
}

# Options that take a value must be followed by one (not by another option).
require_value() {
    if [ -z "$2" ] || [[ "$2" == -* ]]; then
        fail "Option $1 requires a value."
    fi
}

SETNP=""
SETUNLOCK=""
SIFDIR=""
SETBATCH="false"
NBATCH=3
RERUN=false
SETSLURM="false"
PARTITION=""
MEM=""

USED_T=false
USED_J=false

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
    -h|--help)
        usage
        exit 0
        ;;
    -v|--version)
        echo "OrganPipe $VERSION"
        exit 0
        ;;
    -c)
        require_value "$1" "$2"
        shift
        CONFIGFILE=$1
        ;;
    -d)
        require_value "$1" "$2"
        shift
        WORKDIR=$1
        ;;
    -t)
        require_value "$1" "$2"
        shift
        THREADS=$1
        USED_T=true
        ;;
    -j)
        require_value "$1" "$2"
        shift
        THREADS=$1
        USED_J=true
        ;;
    -np)
        SETNP="-np"
        ;;
    -unlock)
        SETUNLOCK="--unlock"
        ;;
    -batch)
        SETBATCH=true
        ;;
    -slurm)
        SETSLURM=true
        ;;
    -partition)
        require_value "$1" "$2"
        shift
        PARTITION=$1
        ;;
    -rerun)
        RERUN=true
        ;;
    -sifdir)
        require_value "$1" "$2"
        shift
        SIFDIR=$1
        ;;
    -nhmmer_db)
        require_value "$1" "$2"
        shift
        NHMMER_DB=$1
        ;;
    -mitos_db)
        require_value "$1" "$2"
        shift
        MITOS_DB=$1
        ;;
    -mem)
        require_value "$1" "$2"
        shift
        MEM=$1
        ;;
    -nbatch)
        require_value "$1" "$2"
        shift
        NBATCH=$1
        ;;
    *)
        fail "Unknown option: $1"
        ;;
    esac
    shift
done

[ -z "$CONFIGFILE" ] && fail "Missing required option -c <configfile>."
[ -z "$WORKDIR" ] && fail "Missing required option -d <workdir>."
if [ "$SETSLURM" != true ] && [ "$USED_T" = false ] && [ "$USED_J" = false ]; then
    fail "Missing required option -t <cores>."
fi

if [ "$SETSLURM" = true ]; then
    [ "$USED_T" = true ] && fail "-t is not allowed with -slurm. Use -j instead."
    [ "$USED_J" = false ] && fail "When using -slurm, you must use -j (not -t) to set the number of jobs submitted to the queue."
    [ -z "$PARTITION" ] && fail "-partition is required when using -slurm."
    [ -n "$MEM" ] && fail "-mem is not allowed with -slurm. Set the memory of each rule in config/resources.yaml."
fi

print_banner

CONFIGFILE=$(realpath "$CONFIGFILE")
SCRIPTDIR="$(dirname "$(readlink -f "$0")")"
WORKDIR=$(realpath "$WORKDIR")

# Paths given on the command line are relative to where OrganPipe.sh was called from
[ -n "$SIFDIR" ] && SIFDIR=$(realpath "$SIFDIR")
[ -n "$NHMMER_DB" ] && NHMMER_DB=$(realpath "$NHMMER_DB")
[ -n "$MITOS_DB" ] && MITOS_DB=$(realpath "$MITOS_DB")

# The scripts and the rules use paths relative to the working directory
[ -d "$WORKDIR" ] || fail "Working directory not found: $WORKDIR"
cd "$WORKDIR"

# MITOS2 runs in a container that only sees the working directory
if [ -n "$MITOS_DB" ]; then
    [[ "$MITOS_DB" == "$WORKDIR"/* ]] || fail "-mitos_db must be inside the working directory ($WORKDIR)."
    MITOS_DB=${MITOS_DB#"$WORKDIR"/}
else
    MITOS_DB="resources/refseq89m"
fi
[ -d "$MITOS_DB" ] || fail "MITOS2 reference directory not found: $MITOS_DB"

if [[ "$CONFIGFILE" == *.csv ]]; then
    CONFIGTYPE="csv"
elif [[ "$CONFIGFILE" == *.yaml ]]; then
    CONFIGTYPE="yaml"
else
    echo "Error: The file is not a valid configuration file (must be .csv or .yaml): $CONFIGFILE"
    exit 1
fi

echo "CONFIGFILE: $CONFIGFILE"
echo "CONFIGTYPE: $CONFIGTYPE"

python $WORKDIR/workflow/scripts/check_config.py -c $CONFIGFILE

if [ "$RERUN" = true ]; then
    python $WORKDIR/workflow/scripts/rerun.py
fi

python "$SCRIPTDIR"/workflow/scripts/create_snakemake_config.py \
        --configfile "$CONFIGFILE"

if [ -n "$SIFDIR" ]; then
    python $WORKDIR/workflow/scripts/singularity.py --sifdir $SIFDIR
    grep -qxF "sif_dir: '$SIFDIR'" config/snakemake_config.yaml || echo "sif_dir: '$SIFDIR'" >> config/snakemake_config.yaml

else
    python $WORKDIR/workflow/scripts/singularity.py --sifdir $WORKDIR/resources/sif_dir
    grep -qxF "sif_dir: 'resources/sif_dir'" config/snakemake_config.yaml || echo "sif_dir: 'resources/sif_dir'" >> config/snakemake_config.yaml
fi

if [ -n "$NHMMER_DB" ]; then
    grep -qxF "nhmmer_db: '$NHMMER_DB'" config/snakemake_config.yaml || echo "nhmmer_db: '$NHMMER_DB'" >> config/snakemake_config.yaml

else
    grep -qxF "nhmmer_db: 'resources/rfam.hmm'" config/snakemake_config.yaml || echo "nhmmer_db: 'resources/rfam.hmm'" >> config/snakemake_config.yaml
fi

grep -qxF "mitos_db: '$MITOS_DB'" config/snakemake_config.yaml || echo "mitos_db: '$MITOS_DB'" >> config/snakemake_config.yaml

mkdir -p $WORKDIR/tmp $WORKDIR/singularity

export SINGULARITY_CACHEDIR=$WORKDIR/singularity
export SINGULARITY_TMPDIR=$WORKDIR/tmp
export TMPDIR=$WORKDIR/tmp

# Build the Snakemake profile from the templates in config/
if [ "$SETSLURM" = true ]; then
    PROFILE=$WORKDIR/profiles/slurm
    CORES=""
    cat $WORKDIR/config/slurm_params.yaml $WORKDIR/config/threads.yaml $WORKDIR/config/resources.yaml | \
        sed "s|{WORKDIR}|$WORKDIR|g; s|{THREADS}|$THREADS|g; s|{PARTITION}|$PARTITION|g" > $PROFILE/config.yaml
else
    PROFILE=$WORKDIR/profiles/local
    CORES="--cores $THREADS --use-singularity"
    [ -n "$MEM" ] && CORES="$CORES --resources mem_mb=$((MEM * 1000))"
    cat $WORKDIR/config/local_params.yaml $WORKDIR/config/threads.yaml | \
        sed "s|{WORKDIR}|$WORKDIR|g; s|{THREADS}|$THREADS|g" > $PROFILE/config.yaml
fi

if [ "$SETBATCH" = true ]; then
    for ((batch=1; batch<=NBATCH; batch++))
    do
        snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile $CORES \
            --scheduler greedy --profile $PROFILE/ \
            $SETNP $SETUNLOCK --batch all=$batch/$NBATCH
    done

else
    snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile $CORES \
        --scheduler greedy --profile $PROFILE/ $SETNP $SETUNLOCK
fi
