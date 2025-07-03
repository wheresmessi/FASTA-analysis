from Bio import SeqIO
from Bio.Seq import Seq
import os

def verify_cds(seq_record):
    seq = str(seq_record.seq).upper()
    errors = []

    if len(seq) <= 300:
        errors.append(f"Length {len(seq)} is not >300 bp")
    if not seq.startswith("ATG"):
        errors.append("Does not start with ATG")
    if seq[-3:] not in ["TAA", "TAG", "TGA"]:
        errors.append("Does not end with valid stop codon")
    if len(seq) % 3 != 0:
        errors.append("Length is not divisible by 3")

    translated_seq = str(Seq(seq).translate(to_stop=False))
    if "*" in translated_seq[:-1]:
        errors.append("Contains internal stop codon")
    
    valid = (len(errors) == 0)
    return valid, errors

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
