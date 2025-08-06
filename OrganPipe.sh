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
SETSLURM="false"

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
    -slurm)
        SETSLURM=true
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

if [ "$SETBATCH" = true ] && [ "$SETSLURM" = true ]; then
    export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
    export TMPDIR=$WORKDIR/tmp && \
    sed -i "s|<WORKDIR>|$WORKDIR|g" $WORKDIR/profiles/slurm/config.yaml && \
    for ((batch=1; batch<=NBATCH; batch++))
    do
        snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile --cores $THREADS \
            --scheduler greedy --rerun-incomplete --profile $WORKDIR/profiles/slurm/ $SETNP $SETUNLOCK --batch all=$batch/$NBATCH
    done

elif [ "$SETBATCH" = true ]; then
    for ((batch=1; batch<=NBATCH; batch++))
    do
        export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
        export TMPDIR=$WORKDIR/tmp && \
        snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile --cores $THREADS --use-singularity \
            --singularity-args "-B $WORKDIR:/mnt -B $WORKDIR/tmp:/tmp --pwd /mnt --no-home --writable" \
            --scheduler greedy --rerun-incomplete $SETNP $SETUNLOCK --batch all=$batch/$NBATCH
    done

elif [ "$SETSLURM" = true ]; then
    export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
    export TMPDIR=$WORKDIR/tmp && \
    sed -i "s|<WORKDIR>|$WORKDIR|g" $WORKDIR/profiles/slurm/config.yaml && \
    snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile --cores $THREADS \
        --scheduler greedy --rerun-incomplete --profile $WORKDIR/profiles/slurm/ $SETNP $SETUNLOCK

else
    export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
    export TMPDIR=$WORKDIR/tmp && \
    snakemake -d $WORKDIR -s $WORKDIR/workflow/Snakefile --cores $THREADS --use-singularity \
        --singularity-args "-B $WORKDIR:/mnt -B $WORKDIR/tmp:/tmp --pwd /mnt --no-home --writable" \
        --scheduler greedy --rerun-incomplete $SETNP $SETUNLOCK
fi
