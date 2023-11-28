import os
import subprocess
import sys

if __name__ == "__main__":
    # Run the command to get the conda bin path
    conda_bin_path = subprocess.check_output("conda info --envs | grep '^GenomeDelta ' | awk '{print $NF}'", shell=True, text=True).strip()

    # Construct the absolute path to "main.sh"
    script_path = os.path.join(conda_bin_path, "bin", "main.sh")

    # Give execute permission and run the script with arguments
    subprocess.call(["chmod", "+x", script_path])

    # Pass command-line arguments to main.sh
    subprocess.call([script_path] + sys.argv[1:])