import argparse
import os
import subprocess

# Parse command-line arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("cluster1", help="")
parser.add_argument("max_distance", help="", type=int)
parser.add_argument("merged_cluster", help="")
args = parser.parse_args()

def get_fasta_files(folder_path):
    files = []
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith(('.fasta', '.fa'))]
    for fa in fasta_files:
        fasta_path = folder_path+"/"+fa
        files.append(fasta_path)
    return files

def extract_info(line):
    info = []
    info.append(line.split(":")[0])
    info.append(int(line.split(":")[1].split("-")[0]))
    info.append(int(line.split(":")[1].split("-")[1]))
    return info

def find_matches(cluster1, cluster2):
    sequences=0
    matches=0
    coupled_list=[]
    with open(cluster1, 'r') as cluster1, open(cluster2, 'r') as cluster2:
        for lineA in cluster1:
            if lineA[0] == ">":
                info = extract_info(lineA)
                chrA = info[0]
                startA = info[1]
                endA = info[2]
                sequences+=1
                for lineB in cluster2:
                    if lineB[0] == ">":
                        info = extract_info(lineB)
                        chrB = info[0]
                        startB = info[1]
                        endB = info[2]
                        if chrA==chrB:
                            if (abs(endA-startB) < args.max_distance) or (abs(startA-endB) < args.max_distance):
                                coupled = (lineA, lineB)
                                coupled_list.append(coupled)
                                matches+=1
                cluster2.seek(0)
    matches_frequency = (matches/sequences)*100
    return matches_frequency, coupled_list

def merge_clusters(coupled, merged):
    with open(merged, 'w') as merged:
        for couple in coupled:
            new_entry = ""
            infoA = extract_info(couple[0])
            infoB = extract_info(couple[1])
            new_start = min(infoA[1], infoB[1])
            new_end = max(infoA[2], infoB[2])
            new_entry = new_entry+infoA[0]+":"+str(new_start)+"-"+str(new_end)+"\n"
            merged.write(new_entry)


cluster_folder = os.path.dirname(args.cluster1)
fasta_files = get_fasta_files(cluster_folder)
checking_cluster = args.cluster1

index = 0
checked_clusters = [args.cluster1]
while index < len(fasta_files):
    cluster = fasta_files[index]
    if cluster not in checked_clusters:
        matches_freq = find_matches(checking_cluster, cluster)
        #print(matches_freq[0])
        #print(matches_freq[1])
        if matches_freq[0] > 50: # CRUCIAL PARAMETER: the minimum percentage of "paired" sequences between the cluster to merge them
            merge_clusters(matches_freq[1], args.merged_cluster)
            checking_cluster = args.merged_cluster
            checked_clusters.append(cluster)
            index = 0  # Restart the loop
        else:
            index += 1  # Move to the next cluster if no merge is performed
    else:
        index += 1  # Move to the next cluster if the current cluster is in the already merged clusters

name_bed = args.merged_cluster
name_bed = name_bed.split(".")[0]
if len(checked_clusters) > 1:
    for cluster in checked_clusters:
        number = cluster.split("/")[-1].split("_")[1].split(".")[0]
        name_bed = name_bed + "_" + str(number)

    bed_folder = os.path.dirname(args.merged_cluster)
    bed_path = name_bed + ".txt"

    # Construct the 'mv' command
    mv_command = ['mv', args.merged_cluster, bed_path]

    # Use subprocess to execute the 'mv' command
    subprocess.run(mv_command, check=True)