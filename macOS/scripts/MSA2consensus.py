import argparse
import os
import statistics

# Parse command-line arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("MSA", help="Path to MSA file")
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

def defragment_MSA(MSA, output):
    with open(MSA, 'r') as input, open(output, 'w') as output:
        i = 0
        for line in input:
            if (line[0] == ">") and (i==0):
                output.write(line)
                i += 1
            elif line[0] == ">":
                output.write("\n"+line)
            else:
                output.write(line[0:-1])
                  

defragment_MSA(args.MSA, args.output+".def")

with open(args.output+".def", 'r') as input, open(args.output, 'w') as output:
    input_basename = os.path.basename(args.output)
    cred = []
    n = 0
    countminus = 0
    scredibility = ""
    for line in input:
        if line.startswith('>'):
            for letter in line:
                if letter=="_":
                    countminus+=1
                elif countminus == 1:
                    scredibility = scredibility + str(letter)
            credibility = float(scredibility)
            countminus = 0
            scredibility = ""
            #print(credibility)
            cred.append(credibility)
            n+=1
    #print(cred)
    cluster_credibility = round(statistics.mean(cred), 2)
    #print(cluster_credibility)
    output.write(">" + input_basename + "_" + str(cluster_credibility) + "_" + str(n) + "\n")
    input.seek(0)
    next(input)
    sequence_length = len(next(input).strip())
    input.seek(0)
    for base in range(0, sequence_length):
        input.seek(0)
        count = [0,0,0,0,0]
        letter = ""
        # A T G C -
        for line in input:
            if line[0] == ">":
                continue
            elif line[base] == "A":
                count[0]+=1
            elif line[base] == "T":
                count[1]+=1
            elif line[base] == "G":
                count[2]+=1
            elif line[base] == "C":
                count[3]+=1
            elif line[base] == "-":
                count[4]+=1
        if max(count) == count[0]:
            letter = "A"
            output.write("A")
        elif max(count) == count[1]:
            letter = "T"
            output.write("T")
        elif max(count) == count[2]:
            letter = "G"
            output.write("G")
        elif max(count) == count[3]:
            letter = "C"
            output.write("C")
        else:
            if count[4]<(count[0]+count[1]+count[2]+count[3]):
                b = [count[0],count[1],count[2],count[3]]
                if max(b) == count[0]:
                    letter = "A"
                    output.write("A")
                elif max(b) == count[1]:
                    letter = "T"
                    output.write("T")
                elif max(b) == count[2]:
                    letter = "G"
                    output.write("G")
                elif max(b) == count[3]:
                    letter = "C"
                    output.write("C")
            else:
                letter = "-"
    output.write("\n")

os.remove(args.output+".def")       