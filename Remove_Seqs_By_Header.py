
"""
Usage:
python Remove_Seqs_By_Header.py input.fasta headers.txt > output.fasta

Obtenido de : https://bioinformatics.stackexchange.com/questions/3931/remove-delete-sequences-by-id-from-multifasta

Misma función pero con awk: 
awk 'BEGIN{while((getline<"ids.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' seq.fa
"""

from Bio import SeqIO
import sys

ffile = SeqIO.parse(sys.argv[1], "fasta")
header_set = set(line.strip() for line in open(sys.argv[2]))

for seq_record in ffile:
    try:
        header_set.remove(seq_record.name)
    except KeyError:
        print(seq_record.format("fasta"))
        continue
if len(header_set) != 0:
    print(len(header_set),'of the headers from list were not identified in the input fasta file.', file=sys.stderr)
