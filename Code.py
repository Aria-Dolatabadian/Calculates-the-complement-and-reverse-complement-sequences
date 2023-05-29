from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
# Define the length of the sequence
seq_length = 50
# Define the nucleotides to choose from
nucleotides = ["A", "C", "G", "T"]
# Generate the random sequence
seq = "".join([random.choice(nucleotides) for _ in range(seq_length)])
# Print the sequence
print(seq)

seq_record = SeqRecord(Seq(seq))

# Add metadata to the sequence record
seq_record.id = "Sequence"
seq_record.description = "A random DNA sequence"

# Write the sequence record to a FASTA file
with open("Sequence.fasta", "w") as output_handle:
    SeqIO.write(seq_record, output_handle, "fasta")


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def complement_sequence(sequence):
    complement = sequence.complement()
    return complement

def reverse_complement_sequence(sequence):
    reverse_complement = sequence.reverse_complement()
    return reverse_complement

def export_fasta(sequence, filename):
    record = SeqRecord(sequence, id='sequence', description='')
    with open(filename, 'w') as file:
        SeqIO.write(record, file, 'fasta')

# Read input FASTA file
input_file = 'Sequence.fasta'
sequences = SeqIO.parse(input_file, 'fasta')

for sequence in sequences:
    seq_id = sequence.id
    seq = sequence.seq

    # Calculate complement sequence
    complement_seq = complement_sequence(seq)

    # Calculate reverse complement sequence
    reverse_complement_seq = reverse_complement_sequence(seq)

    # Export complement sequence to FASTA file
    complement_filename = f'{seq_id}_complement.fasta'
    export_fasta(complement_seq, complement_filename)

    # Export reverse complement sequence to FASTA file
    reverse_complement_filename = f'{seq_id}_reverse_complement.fasta'
    export_fasta(reverse_complement_seq, reverse_complement_filename)

    print(f'Complement sequence exported to: {complement_filename}')
    print(f'Reverse complement sequence exported to: {reverse_complement_filename}')
