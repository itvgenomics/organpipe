from Bio import SeqIO
import os


def extract_sequences(genbank_file, fasta_file, featuretype):
    records = SeqIO.parse(genbank_file, "genbank")
    cds_sequences = []

    for record in records:
        for feature in record.features:
            if feature.type == featuretype:
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                nucleotide_seq = feature.location.extract(record).seq
                cds_sequences.append((product, nucleotide_seq))

    with open(fasta_file, "w") as output:
        for i, (product, sequence) in enumerate(cds_sequences, start=1):
            header = f">{str(product).replace(' ', '_')}"
            output.write(f"{header}\n{sequence}\n")


def split_fasta(fasta_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    index = 1
    # Iterate over each record in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        output_path = os.path.join(output_dir, f"{index}-{record.id}.fasta")

        # Write each record to a separate file
        with open(output_path, "w") as outfile:
            SeqIO.write(record, outfile, "fasta")

        index += 1


def get_fasta_headers(fasta_file):
    headers = []
    index = 1
    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(str(f"{index}-{record.id}"))
        index += 1
    return headers


def extract_fasta_headers(fasta_file, output_file):
    headers = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.id)
    with open(output_file, "w") as f:
        for header in headers:
            f.write(header + "\n")
    return headers


# Function to read headers from a file
def read_headers(header_file):
    with open(header_file, "r") as f:
        return [line.strip() for line in f]


def main(featuretype, genbank_file, output_file):
    # output_file = f"{str(genbank_file).rsplit('.', 1)[0]}_{featuretype}.fasta"
    extract_sequences(genbank_file, output_file, featuretype)
    return output_file
