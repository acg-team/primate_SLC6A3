from Bio import SeqIO

def main():
    fasta_file = "data/primates_SLC6A3.fa"
    sequences = SeqIO.parse(fasta_file, "fasta")

    for record in SeqIO.parse(fasta_file, "fasta"):        
        print(record.id, len(record))
    
    seqs_renamed = [record for record in SeqIO.parse(fasta_file, "fasta")]
    for seq in seqs_renamed:
        seq.id = seq.id.replace("/", "-")
        seq.description = seq.description.replace("/", "-")
    for s in seqs_renamed:
        print(len(s))
    
    SeqIO.write(seqs_renamed, "data/primates_SLC6A3_reformat.fa", "fasta")

if __name__ == "__main__":
    main()
