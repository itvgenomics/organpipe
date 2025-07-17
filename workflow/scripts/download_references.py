#!/usr/bin/env python

"""
This script finds mitogenomes from related species in public databases.

License:
    Copyright 2022 Ksenia Krasheninnikova
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

from argparse import ArgumentParser
import os
from Bio import Entrez
from Bio import SeqIO
from io import StringIO
import logging
import sys

FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

logging.basicConfig(
    level=logging.INFO,
    stream=sys.stdout,
    format=FORMAT,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def get_lineage(species):
    handle = Entrez.esearch(db="taxonomy", term=species, idtype="acc")
    record = Entrez.read(handle)
    if not len(record["IdList"]):
        raise Exception("No such species in NCBI!")
    if len(record["IdList"]) > 1:
        raise Exception("More than one appropriate entries in NCBI!")
    handle = Entrez.efetch(db="Taxonomy", id=record["IdList"][0], retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    for e in records[0]["Lineage"].split(" ")[::-1]:
        yield e


def find_full_mito(
    group, outfolder, length_threshold, considered, org_type="mitochondrion", n=1
):
    term = (
        '("'
        + group
        + '"[Organism] AND complete '
        + "genome[All Fields]) AND "
        + org_type
        + "[filter]]"
    )
    handle = Entrez.esearch(db="nucleotide", term=term, idtype="acc")
    record = Entrez.read(handle)
    if record["IdList"]:
        for ncbi_code in record["IdList"]:
            handle = Entrez.efetch(
                db="nucleotide", id=ncbi_code, rettype="gb", retmode="text"
            )
            record_ = handle.read()
            handle.close()
            for seqrecord in SeqIO.parse(StringIO(record_), "gb"):
                if seqrecord.id in considered:
                    continue
                considered.add(seqrecord.id)
                if len(seqrecord) > length_threshold:
                    if len(seqrecord.features) < 10:
                        logging.info(
                            "Not enough features in gb file! skipping " + seqrecord.id
                        )
                        continue
                    with open(os.path.join(outfolder, ncbi_code + ".gb"), "w") as out:
                        out.write(record_)
                    handle = Entrez.efetch(
                        db="nucleotide", id=ncbi_code, rettype="fasta", retmode="text"
                    )
                    record_ = handle.read()
                    handle.close()
                    with open(
                        os.path.join(outfolder, ncbi_code + ".fasta"), "w"
                    ) as out:
                        out.write(record_)
                    logging.info(
                        "output is written to "
                        + os.path.join(outfolder, ncbi_code)
                        + ".[gb,fasta]"
                    )
                    n -= 1
                    if n == 0:
                        return considered, n
    return considered, n


def main(species, email, outfolder, min_length, n):
    org_type = "mitochondrion"
    if n < 1:
        logging.info("Number of genomes to report must be at least 1 (default)")
        return

    Entrez.email = email
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder, exist_ok=True)

    logging.info("Looking for " + org_type + " for " + species)
    considered = set()
    for g in [species] + list(get_lineage(species)):
        if n > 0:
            logging.info("Looking for an appropriate organelle among " + g)
            considered, n = find_full_mito(
                g, outfolder, min_length, considered, org_type, n
            )

    if n == 1:  # Assuming n was originally set to 1
        logging.info("No appropriate mitogenome found")
