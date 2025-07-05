from Bio import SeqIO
from collections import Counter
import os

RESULTS_FOLDER = "results"
INPUT_FILE = "uploads/input_for_analysis.txt"

def calculate_gc(seq):
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq) * 100) if len(seq) > 0 else 0

def calculate_gc_by_codon_position(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
    if not codons:
        return (0, 0, 0)
    gc1 = sum(1 for c in codons if c[0] in "GC")
    gc2 = sum(1 for c in codons if c[1] in "GC")
    gc3 = sum(1 for c in codons if c[2] in "GC")
    total = len(codons)
    return (gc1 / total * 100, gc2 / total * 100, gc3 / total * 100)

def run_combined_analysis():
    records = list(SeqIO.parse(INPUT_FILE, "fasta"))
    base_summary = []
    csv_summary = []

    header = ("Species\tNum_CDS\tA_count\tT_count\tC_count\tG_count\tA%\tT%\tC%\tG%\t"
              "A3_count\tT3_count\tC3_count\tG3_count\tA3%\tT3%\tC3%\tG3%\tL_aa\t"
              "Overall_GC%\tGC1%\tGC2%\tGC3%")
    csv_header = header.replace("\t", ",")

    base_summary.append(header)
    csv_summary.append(csv_header)

    species = os.path.basename(INPUT_FILE).replace("_valid.fasta", "").replace(".txt", "")

    overall_counts = Counter()
    third_counts = Counter()
    total_aa = 0
    gc_total, gc1_list, gc2_list, gc3_list = [], [], [], []

    for r in records:
        seq = str(r.seq).upper()
        overall_counts.update(seq)

        codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
        for codon in codons:
            third_counts[codon[2]] += 1

        aa_len = (len(seq) // 3) - 1 if seq[-3:] in ["TAA", "TAG", "TGA"] else len(seq) // 3
        total_aa += aa_len

        gc_total.append(calculate_gc(seq))
        gc1, gc2, gc3 = calculate_gc_by_codon_position(seq)
        gc1_list.append(gc1)
        gc2_list.append(gc2)
        gc3_list.append(gc3)

    total_bases = sum(overall_counts[b] for b in "ATGC")
    total_third = sum(third_counts[b] for b in "ATGC")

    def pct(n): return (n / total_bases * 100) if total_bases else 0
    def pct3(n): return (n / total_third * 100) if total_third else 0

    row = (
        f"{species}\t{len(records)}\t"
        f"{overall_counts['A']}\t{overall_counts['T']}\t{overall_counts['C']}\t{overall_counts['G']}\t"
        f"{pct(overall_counts['A']):.2f}\t{pct(overall_counts['T']):.2f}\t{pct(overall_counts['C']):.2f}\t{pct(overall_counts['G']):.2f}\t"
        f"{third_counts['A']}\t{third_counts['T']}\t{third_counts['C']}\t{third_counts['G']}\t"
        f"{pct3(third_counts['A']):.2f}\t{pct3(third_counts['T']):.2f}\t{pct3(third_counts['C']):.2f}\t{pct3(third_counts['G']):.2f}\t"
        f"{total_aa}\t"
        f"{sum(gc_total)/len(gc_total):.2f}\t"
        f"{sum(gc1_list)/len(gc1_list):.2f}\t"
        f"{sum(gc2_list)/len(gc2_list):.2f}\t"
        f"{sum(gc3_list)/len(gc3_list):.2f}"
    )

    base_summary.append(row)
    csv_summary.append(row.replace("\t", ","))

    with open(os.path.join(RESULTS_FOLDER, "combined_summary.txt"), "w") as txtf:
        txtf.write("\n".join(base_summary))

    with open(os.path.join(RESULTS_FOLDER, "combined_summary.csv"), "w") as csvf:
        csvf.write("\n".join(csv_summary))