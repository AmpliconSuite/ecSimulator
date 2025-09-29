"""
File for basic utilities.
"""
from math import ceil


def get_fasta_length(fasta):
    tlen = 0
    with open(fasta) as infile:
        for line in infile:
            if line.startswith(">"):
                continue

            tlen+=len(line.rstrip())

    print("Fasta size is " + str(tlen) + "bp")
    return tlen


def pseudocircularize_fasta(fasta, coverage):
    # read a fasta and makes c copies of each entry concatenated with itself
    # writes a new fasta and return the path of that file
    base, ext = os.path.splitext(fasta)
    output_fasta = f"{base}_circularized{ext}"
    c = max(2, int(ceil(coverage)))
    with open(fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        seq_name = ""
        seq_data = ""

        for line in infile:
            if line.startswith('>'):
                # Write the previous entry if it exists
                if seq_name and seq_data:
                    pseudocircularized_seq = seq_data * c
                    outfile.write(f"{seq_name}\n")
                    outfile.write(f"{pseudocircularized_seq}\n")

                # Start a new entry
                seq_name = line.strip()
                seq_data = ""
            else:
                seq_data += line.strip()

        # Write the last entry
        if seq_name and seq_data:
            pseudocircularized_seq = seq_data * c
            outfile.write(f"{seq_name}\n")
            outfile.write(f"{pseudocircularized_seq}\n")

    return output_fasta


def compute_number_of_reads_to_simulate(coverage, fasta_len, mean_rl):
    nreads = int(ceil(coverage * fasta_len / mean_rl))
    print("Number of reads to simulate for {}x coverage: {}".format(str(coverage), str(nreads)))
    return nreads
