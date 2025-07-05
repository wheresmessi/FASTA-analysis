from Bio import SeqIO
from Bio.Seq import Seq

def verify_cds(seq_record):
    seq = str(seq_record.seq).upper()
    errors = []

    if len(seq) <= 300:
        errors.append("Length <= 300 bp")
    if not seq.startswith("ATG"):
        errors.append("Does not start with ATG")
    if seq[-3:] not in ["TAA", "TAG", "TGA"]:
        errors.append("Invalid stop codon")
    if len(seq) % 3 != 0:
        errors.append("Not divisible by 3")

    translated = str(Seq(seq).translate(to_stop=False))
    if "*" in translated[:-1]:
        errors.append("Contains internal stop codon")

    return (len(errors) == 0), errors

def validate_fasta_file(filepath):
    valid_records = []
    invalids = []

    for record in SeqIO.parse(filepath, "fasta"):
        valid, errors = verify_cds(record)
        if valid:
            valid_records.append(record)
        else:
            invalids.append({"id": record.id, "errors": errors})

    return valid_records, invalids