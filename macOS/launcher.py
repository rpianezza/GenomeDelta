import os
import subprocess
import sys
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="GenomeDelta launcher")
parser.add_argument("--main", help="Path to main file")
args = parser.parse_args()

if __name__ == "__main__":
    script_path = args.main

    # Give execute permission and run the script with arguments
    subprocess.call(["chmod", "+x", script_path])

    # Pass command-line arguments to main.sh
    subprocess.call([script_path] + sys.argv[1:])