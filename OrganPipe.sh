#!/bin/bash
set -e

cat << 'EOF'

### OrganPipe: An automated tool to facilitate the assembly, annotation, and curation of mitochondrial and chloroplast genomes
#Authors: Renato R. Moreira-Oliveira, Bruno Marques Silva, Michele Molina, Marx Oliveira-Lima, Tiago Ferreira Leão, Santelmo Vasconcelos, Gisele Lopes Nunes. 2024

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

SETNP=""
SETUNLOCK=""
SIFDIR=""
SETBATCH="false"
NBATCH=15
RERUN=false
SETNOTEMP=""

while [ "$1" != "" ]; do
    case $1 in
    -c)
        shift
        CONFIGFILE=$1
        ;;
    -d)
        shift
        WORKDIR=$1
        ;;
    -t)
        shift
        THREADS=$1
        ;;
    -np)
        SETNP="-np"
        ;;
    -unlock)
        SETUNLOCK="--unlock"
        ;;
    -notemp)
        SETNOTEMP="--notemp"
        ;;
    -batch)
        SETBATCH=true
        ;;
    -rerun)
        RERUN=true
        ;;
    -sifdir)
        shift
        SIFDIR=$1
        ;;
    -nbatch)
        shift
        NBATCH=$1
        ;;
    *)
        echo "Invalid option: $1"
        echo "Usage: bash OrganPipe.sh -d </path/to/work/dir> -t <n_threads> -c </path/to/configfile> [options]"
        echo ""
        echo "Required:"
        echo "  -d </path/to/work/dir>    Path to your working directory where all the workflow files are"
        echo "  -c </path/to/config.yaml> Overwrite the default configuration file (e.g. config/config.yaml)"
        echo "  -t {int}                  Number of threads to use"
        echo ""
        echo "Optional Flags:"
        echo "  -np                       Perform a dry run to preview job execution without running"
        echo "  -unlock                   Unlock the working directory if Snakemake is locked"
        echo "  -batch                    Improve DAG resolution time for large workflows"
        echo "  -nbatch {int}             Set batch size (default: 15) for large sample sets"
        echo "  -sifdir </path>           Directory to build/store Singularity images (default: resources/sif_dir)"
        echo "  -rerun                    Clean previous results/temp files for clean reprocessing"
        echo "  -notemp                   Prevent deletion of temporary files (e.g. for partial runs)"
        exit 1
        ;;
    esac
    shift
done


CONFIGFILE=$(realpath "$CONFIGFILE")
SCRIPTDIR="$(dirname "$(readlink -f "$0")")"
WORKDIR=$(realpath "$WORKDIR")

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
    SIFDIR=$(realpath "$SIFDIR")
    python $WORKDIR/workflow/scripts/singularity.py --sifdir $SIFDIR
    grep -qxF "sif_dir: '$SIFDIR'" config/snakemake_config.yaml || echo "sif_dir: '$SIFDIR'" >> config/snakemake_config.yaml

else
    python $WORKDIR/workflow/scripts/singularity.py --sifdir $WORKDIR/resources/sif_dir
    grep -qxF "sif_dir: 'resources/sif_dir'" config/snakemake_config.yaml || echo "sif_dir: 'resources/sif_dir'" >> config/snakemake_config.yaml
fi

mkdir -p $WORKDIR/tmp $WORKDIR/singularity

if [ "$SETBATCH" = true ]; then
    for ((batch=1; batch<=NBATCH; batch++))
    do
        export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
        export TMPDIR=$WORKDIR/tmp && \
        snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile --cores $THREADS --use-singularity \
            --singularity-args "-B $WORKDIR:/mnt -B $WORKDIR/tmp:/tmp --pwd /mnt --no-home --writable" \
            --scheduler greedy --rerun-incomplete $SETNP $SETUNLOCK $SETNOTEMP --batch all=$batch/$NBATCH
    done

else
    export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
    export TMPDIR=$WORKDIR/tmp && \
    snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile --cores $THREADS --use-singularity \
        --singularity-args "-B $WORKDIR:/mnt -B $WORKDIR/tmp:/tmp --pwd /mnt --no-home --writable" \
        --scheduler greedy --rerun-incomplete $SETNP $SETUNLOCK $SETNOTEMP
fi
