from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.Seq import Seq
import os
import logging
import sys
import re
import subprocess
from naming import set_genbank_identity, warn, write_genbank

# Sequence ID given to table2asn. The LOCUS name is replaced afterwards, so its
# length limit does not apply to the names written by OrganPipe.
TABLE2ASN_ID = "Asm_Contig"

# CPGAVAS2 results are named after its job id: {pid}.gbf, {pid}.html...
CPGAVAS2_GENBANK = re.compile(r"(\d+)\.gbf")

GENE_FEATURES = {"gene", "CDS", "tRNA", "rRNA"}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create GenBank file from the circularized mtGenome"
    )
    parser.add_argument(
        "--software", help="Software used (chloe or cpgavas2)", required=True,
    )
    parser.add_argument(
        "--annotation_dir", help="Annotation directory of the assembly", required=True,
    )
    parser.add_argument("--output", help="GenBank file to write", required=True)
    parser.add_argument(
        "--gene2product",
        type=str,
        nargs="?"
    )
    parser.add_argument("--fasta", help="Annotated assembly (Chloe)")
    parser.add_argument("--name", help="LOCUS name of the GenBank file")
    parser.add_argument("--definition", help="DEFINITION of the GenBank file")

    args = parser.parse_args()

    return args

def create_features_table(gene2product, gff_file, features_file, fasta):
    D_gene2product = {}
    for line in open(gene2product):
        line = line.strip()
        if line and line[0] != "#":
            D_gene2product[line.split("\t")[0]] = line.split("\t")[1]

    GFF_File = gff_file
    outFiveColumnFile = features_file
    genomeFile = fasta

    record = SeqIO.read(genomeFile, "fasta")
    seqID = record.id

    L = open(GFF_File).readlines()

    contig_line = ""

    annotate_block_start = 0
    annotate_block_end = 0

    for i in range(len(L)):
        line = L[i]
        if "##sequence-region" in line:
            contig_line = line
            annotate_block_start = i + 1

        if line == "##FASTA\n" and annotate_block_end == 0:
            annotate_block_end = i - 1
    annotate_block = L[annotate_block_start : annotate_block_end + 1]

    contigLength = int(contig_line.split("\t")[3])
    (D_gene, D_mRNA) = GFF_Parse_GetFeatureBlock(annotate_block)

    # outFiveColumnFile = sys.argv[2]
    # seqID = sys.argv[3]
    out = open(outFiveColumnFile, "w")
    out.write(">Feature Asm_Contig\n")

    for key in D_gene.keys():
        geneID = key
        genename = D_gene[key]["genename"]

        # geneID = key
        # genename = D_gene[key]["genename"]
        geneLine = D_gene[key]["geneLine"]
        if genename not in D_gene2product.keys():
            print("Warning:\t", genename, "not in D_gene2product")
            product = "Unknown"
        else:
            product = D_gene2product[genename]

        (gene_start, gene_end, gene_strand) = (
            geneLine.split("\t")[3],
            geneLine.split("\t")[4],
            geneLine.split("\t")[6],
        )

        mRNALine = D_gene[key]["mRNALine"]
        mRNAID = D_gene[key]["mRNAID"]
        (mRNA_start, mRNA_end, mRNA_strand) = (
            mRNALine.split("\t")[3],
            mRNALine.split("\t")[4],
            mRNALine.split("\t")[6],
        )

        featuretype = D_mRNA[mRNAID]["featuretype"]
        L_exon = D_mRNA[mRNAID]["exon"]

        L_CDS = D_mRNA[mRNAID]["CDS"]

        if genename not in ["rps12", "trnH"]:

            if gene_strand == "+":
                out.write(gene_start + "\t" + gene_end + "\tgene\t\t\n")
            else:
                out.write(gene_end + "\t" + gene_start + "\tgene\t\t\n")

            out.write("\t\t\tgene\t" + genename + "\n")

        else:
            # contigLength
            if genename == "trnH":
                if int(gene_start) < int(gene_end):

                    if gene_strand == "+":
                        out.write(gene_start + "\t" + gene_end + "\tgene\t\t\n")
                    else:
                        out.write(gene_end + "\t" + gene_start + "\tgene\t\t\n")
                    out.write("\t\t\tgene\t" + genename + "\n")

                else:
                    # gene_start = int(gene_start)
                    if contigLength != int(gene_start):
                        if contigLength - int(gene_start) < 100:
                            out.write(
                                gene_start + "\t" + str(contigLength) + "\tgene\t\t\n"
                            )
                            out.write("1" + "\t" + gene_end + "\t\t\t\n")
                            out.write("\t\t\tgene\t" + genename + "\n")

                    else:
                        out.write("1" + "\t" + gene_end + "\tgene\t\t\n")
                        out.write("\t\t\tgene\t" + genename + "\n")

        if featuretype == "CDS":

            L_CDS_Pos = []

            for CDSLine in L_CDS:
                (CDS_start, CDS_end, CDS_strand) = (
                    CDSLine.split("\t")[3],
                    CDSLine.split("\t")[4],
                    CDSLine.split("\t")[6],
                )
                L_CDS_Pos.append([CDS_start, CDS_end, CDS_strand])

            if genename == "rps12":

                if len(L_CDS_Pos) == 1:

                    if gene_strand == "+":
                        out.write(gene_start + "\t" + gene_end + "\tgene\t\t\n")
                    else:
                        out.write(gene_end + "\t" + gene_start + "\tgene\t\t\n")

                    out.write("\t\t\tgene\t" + genename + "\n")

                else:
                    exon1_pos_part = L_CDS_Pos[0]
                    L_exonOther_pos_part = L_CDS_Pos[1:]

                    (exon1_start, exon1_end, exon1_strand) = exon1_pos_part
                    if exon1_strand == "+":
                        out.write(exon1_start + "\t" + exon1_end + "\tgene\t\t\n")
                    else:
                        out.write(exon1_end + "\t" + exon1_start + "\tgene\t\t\n")

                    exonOther_strand = L_exonOther_pos_part[0][2]

                    if exonOther_strand == "+":
                        exonOther_start = L_exonOther_pos_part[0][0]
                        exonOther_end = L_exonOther_pos_part[-1][1]
                        out.write(exonOther_start + "\t" + exonOther_end + "\t\t\t\n")
                    else:
                        exonOther_start = L_exonOther_pos_part[-1][0]
                        exonOther_end = L_exonOther_pos_part[0][1]
                        out.write(exonOther_end + "\t" + exonOther_start + "\t\t\t\n")

                    out.write("\t\t\tgene\t" + genename + "\n")
                    out.write("\t\t\ttrans_splicing\t\n")
            first_CDS_flag = 1
            for CDS_Pos_Part in L_CDS_Pos:
                (CDS_start, CDS_end, CDS_strand) = CDS_Pos_Part
                if first_CDS_flag:
                    if CDS_strand == "+":
                        # L_CDS_line.append(CDS_start+"\t"+CDS_end+"\tCDS\t\t\n")
                        out.write(CDS_start + "\t" + CDS_end + "\tCDS\t\t\n")

                    else:
                        # L_CDS_line.append(CDS_end+"\t"+CDS_start+"\tCDS\t\t\n")
                        out.write(CDS_end + "\t" + CDS_start + "\tCDS\t\t\n")

                    first_CDS_flag = 0

                else:
                    if CDS_strand == "+":
                        # L_CDS_line.append(CDS_start+"\t"+CDS_end+"\t\t\t\n")
                        out.write(CDS_start + "\t" + CDS_end + "\t\t\t\n")

                    else:
                        # L_CDS_line.append(CDS_end+"\t"+CDS_start+"\t\t\t\n")
                        out.write(CDS_end + "\t" + CDS_start + "\t\t\t\n")

            if genename == "rps12":
                out.write("\t\t\tproduct\t" + product + "\n")
                out.write("\t\t\tcodon_start\t1" + "\n")
                if len(L_CDS_Pos) >= 2:
                    out.write("\t\t\texception\ttrans-splicing" + "\n")
            else:
                out.write("\t\t\tproduct\t" + product + "\n")
                out.write("\t\t\tcodon_start\t1" + "\n")

        elif featuretype in ["tRNA", "rRNA"]:
            # genename
            L_exon_Pos = []

            if genename != "trnH":
                for exonLine in L_exon:
                    (exon_start, exon_end, exon_strand) = (
                        exonLine.split("\t")[3],
                        exonLine.split("\t")[4],
                        exonLine.split("\t")[6],
                    )
                    L_exon_Pos.append([exon_start, exon_end, exon_strand])
            else:
                if int(gene_start) < int(gene_end):
                    for exonLine in L_exon:
                        (exon_start, exon_end, exon_strand) = (
                            exonLine.split("\t")[3],
                            exonLine.split("\t")[4],
                            exonLine.split("\t")[6],
                        )
                        L_exon_Pos.append([exon_start, exon_end, exon_strand])

                else:
                    if contigLength != int(gene_start):
                        if contigLength - int(gene_start) < 100:
                            L_exon_Pos = [
                                (gene_start, str(contigLength), "+"),
                                ("1", gene_end, "+"),
                            ]
                    else:
                        L_exon_Pos = [("1", gene_end, "+")]

                # L_exon_Pos.sort(key=lambda x:int(x[0])

            if L_exon_Pos:
                first_exon_flag = 1
                for exon_Pos_Part in L_exon_Pos:
                    (exon_start, exon_end, exon_strand) = exon_Pos_Part
                    if first_exon_flag:
                        if exon_strand == "+":
                            out.write(
                                exon_start
                                + "\t"
                                + exon_end
                                + "\t%s\t\t\n" % featuretype
                            )
                        else:
                            out.write(
                                exon_end
                                + "\t"
                                + exon_start
                                + "\t%s\t\t\n" % featuretype
                            )
                        first_exon_flag = 0

                    else:
                        if exon_strand == "+":
                            out.write(exon_start + "\t" + exon_end + "\t\t\t\n")
                        else:
                            out.write(exon_end + "\t" + exon_start + "\t\t\t\n")
                out.write("\t\t\tproduct\t" + product + "\n")

    out.close()

    print("GenBank File have successfully be created.")


def GFF_Parse_GetFeatureBlock(annotate_block):

    D_gene = {}
    D_mRNA = {}

    for line in annotate_block:
        # line = line.strip()
        if line:

            featuretype = line.split("\t")[2]
            # print featuretype
            if featuretype == "gene":
                geneID = line.split("\t")[8].split("ID=")[1].split(";")[0].strip()
                # print geneID
                D_gene[geneID] = {}
                D_gene[geneID]["geneLine"] = line
                D_gene[geneID]["geneID"] = geneID

            elif featuretype in ["rRNA", "mRNA", "tRNA", "CDS"]:  # Might produce an error if unexpected feature types appear
                mRNAID = line.split("\t")[8].split("ID=")[1].split(";")[0]
                mRNAParent = (
                    line.split("\t")[8].split("parent=")[1].split(";")[0].strip()
                )
                genename = line.split("\t")[8].split("gene=")[1].split(";")[0]

                D_gene[mRNAParent]["mRNALine"] = line.strip()
                D_gene[mRNAParent]["mRNAID"] = mRNAID
                D_gene[mRNAParent]["genename"] = genename

                if mRNAID not in D_mRNA:
                    D_mRNA[mRNAID] = {}
                    D_mRNA[mRNAID]["exon"] = []
                    D_mRNA[mRNAID]["CDS"] = []

                D_mRNA[mRNAID]["exon"].append(line.strip())
                D_mRNA[mRNAID]["CDS"].append(line.strip())

                D_mRNA[mRNAID]["featuretype"] = featuretype
                D_mRNA[mRNAID]["genename"] = genename

    return (D_gene, D_mRNA)

def write_table2asn_fasta(assembly_fasta, fasta_file):
    """Copy the assembly with the sequence ID that table2asn expects."""
    record = SeqIO.read(assembly_fasta, "fasta")
    record.id = TABLE2ASN_ID
    record.description = ""
    SeqIO.write(record, fasta_file, "fasta")


def has_genes(record):
    return any(feature.type in GENE_FEATURES for feature in record.features)


def clear_template_organism(record):
    """CPGAVAS2 copies the organism and references of its example file (Arabidopsis thaliana)."""
    record.annotations["source"] = "chloroplast Unknown"
    record.annotations["organism"] = "Unknown."
    record.annotations["taxonomy"] = ["Unclassified"]
    record.annotations.pop("references", None)
    for feature in record.features:
        if feature.type == "source":
            feature.qualifiers["organism"] = ["unknown"]
            feature.qualifiers.pop("db_xref", None)


def write_named_genbank(record, output, name, definition):
    if name:
        set_genbank_identity(record, name, definition or record.description)
    write_genbank(record, output)


def chloe_genbank(args):
    """Write the GenBank file of a Chloe annotation. Returns True if it was written."""
    for file in sorted(os.listdir(args.annotation_dir)):
        if not file.endswith(".chloe.gff"):
            continue

        # Chloe runs with --no-transform: the annotation is on the assembly as it is,
        # and Chloe writes no FASTA of it
        base = f"{args.annotation_dir}/{file[:-len('.gff')]}"
        write_table2asn_fasta(args.fasta, f"{base}.fsa")
        create_features_table(fasta=f"{base}.fsa",
                    gff_file=f"{base}.gff",
                    gene2product=args.gene2product,
                    features_file=f"{base}.tbl")

        command = [
            './resources/table2asn.linux64',
            '-j', '[topology=circular] [Completedness=complete] [location=chloroplast] [gcode=11]',
            '-a', 's',
            '-i', f"{base}.fsa",
            '-f', f"{base}.tbl",
            '-V', 'vb'
        ]
        subprocess.run(command, check=True)

        write_named_genbank(SeqIO.read(f"{base}.gbf", "genbank"), args.output, args.name, args.definition)
        return True

    warn(f"Chloe produced no annotation in {args.annotation_dir}, so there is no GenBank file.")
    return False


def cpgavas2_genbank(args):
    """Write the GenBank file of a CPGAVAS2 annotation. Returns True if it was written."""
    for file in sorted(os.listdir(args.annotation_dir)):
        if not CPGAVAS2_GENBANK.fullmatch(file):
            continue

        record = SeqIO.read(f"{args.annotation_dir}/{file}", "genbank")
        if not has_genes(record):
            # CPGAVAS2 still writes a GenBank file, with the source feature only
            warn(
                f"CPGAVAS2 annotated no genes in {args.annotation_dir}, so there is no GenBank file. "
                "See its run_cpgavas2.log."
            )
            return False

        clear_template_organism(record)
        write_named_genbank(record, args.output, args.name, args.definition)
        return True

    warn(f"CPGAVAS2 produced no annotation in {args.annotation_dir}, so there is no GenBank file.")
    return False


if __name__ == "__main__":
    args = parse_arguments()

    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%d-%m-%Y %H:%M:%S",
    )

    logging.info("Running genbank.py")

    logging.info(
        f"""
-   ANNOTATION: {args.annotation_dir}
-   SOFTWARE: {args.software}
"""
    )

    if args.software == "chloe":
        chloe_genbank(args)
    elif args.software == "cpgavas2":
        cpgavas2_genbank(args)
