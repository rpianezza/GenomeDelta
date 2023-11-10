import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Assign credibility scores to each gap.")
parser.add_argument("bed", help="Path to bed file")
parser.add_argument("cred", help="Path to credibility file")
parser.add_argument("fai", help="Path to fai file")
args = parser.parse_args()



def find_flanking(chr, start, end, fai, bases):
    flanking = []
    with open(fai, 'r') as fai_file:
        for line in fai_file:
            if chr == line.split('\t')[0]:
                end_chr = int(line.split('\t')[1])
    if bases < 0:
        start_flank = start+int(bases)
        if start_flank <= 0:
            start_flank = 1
        end_flank = start
    else:
        start_flank = end
        end_flank = end+int(bases)
        if end_flank > end_chr:
            end_flank = end_chr 
    flanking.append(chr)
    flanking.append(str(start_flank))
    flanking.append(str(end_flank))
    return flanking

def extract_credibility(interval, cred):
    credibility = 0
    with open(cred, 'r') as credibility_file:
        for line in credibility_file:
            if (line.split(' ')[0]==interval[0]) and (line.split(' ')[1]==interval[1]) and (line.split(' ')[2]==interval[2]):
                credibility = round(float(line.split(' ')[4]),2)
    #if credibility == 0:
        #print("Flanking interval not found for " + str(interval) + " (maybe is telomeric?)")
    return credibility

output_file_name = args.bed.replace('.bed', '-credibility.bed')

with open(args.bed, 'r') as bed_file, open(output_file_name, 'w') as output_file:
    for line in bed_file:
        chr = line.split('\t')[0]
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2])
        interval_left = find_flanking(chr, start, end, args.fai, -10000)
        interval_right = find_flanking(chr, start, end, args.fai, 10000)
        left_credibility = extract_credibility(interval_left, args.cred)
        right_credibility = extract_credibility(interval_right, args.cred)
        cred = round(((left_credibility + right_credibility) / 2), 2)
        seq_name = chr + ":" + str(start) + "-" + str(end) + "-" + str(cred)

        # Append the 'cred' value as the 4th column
        new_line = f"{line.strip()}\t{seq_name}\t{cred}\n"
        output_file.write(new_line)