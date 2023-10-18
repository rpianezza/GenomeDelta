import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Find clusters in BLAST matrix.")
parser.add_argument("blast", help="Path to blast file")
parser.add_argument("fasta", help="Path to fasta file")
parser.add_argument("output", help="Path to the output folder")
args = parser.parse_args()

'''

'''
def find_clusters(blast):
    i = 0
    query = ""
    subject = ""
    clusters = [[]]
    with open(blast, 'r') as blast_file:
        for line in blast_file:
            if query == "":
                query = line.split("\t")[0]
                subject = line.split("\t")[1]
                clusters[i].append(query)
                #print("New query: " + str(query))
                #print("New cluster: " + str(clusters[i]))
            elif line.split("\t")[0] != query:
                query = line.split("\t")[0]
                #print("New query: " + str(query))
                if not any(query in sublist for sublist in clusters):
                    i += 1
                    clusters.append([query])
                    #print("New cluster: " + str(clusters[i]))
            else:
                subject = line.split("\t")[1]
                if not any(subject in sublist for sublist in clusters):
                    #print("Appending to cluster: " + str(subject))
                    clusters[i].append(subject)
    return clusters

def filter(fasta, selection, out):
    with open(fasta, 'r') as fasta_file, open(out, 'w') as output_file:
        write_line = False
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                gene_name = line[1:]
                if gene_name in selection:
                    write_line = True
                    output_file.write(line + "\n")  # Write the current line to output
            elif write_line:
                output_file.write(line + "\n")  # Write the next line to output
                write_line = False

def check_len(fasta):
    with open(fasta, 'r') as fasta_file:
        lens = []
        saved = []
        seq = ""
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                continue
            else:
                sequence_len = len(line)
                lens.append(sequence_len)
        mean_len = sum(lens)/len(lens)
        fasta_file.seek(0)
        for line in fasta_file:
            if line.startswith(">"):
                seq = line[1:-1]
            else:
                sequence_len = len(line)
                if sequence_len < (mean_len*2):
                    saved.append(seq)
    return(saved)

clusters = find_clusters(args.blast)
for i, cluster in enumerate(clusters):
    if len(cluster) > 2:
        filter(args.fasta, cluster, args.output+"cluster_"+str(i)+".fasta.raw")
        correct_len = check_len(args.output+"cluster_"+str(i)+".fasta.raw")
        filter(args.output+"cluster_"+str(i)+".fasta.raw", correct_len, args.output+"cluster_"+str(i)+".fasta")
        os.remove(args.output+"cluster_"+str(i)+".fasta.raw")