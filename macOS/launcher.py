import os
import subprocess
import sys

if __name__ == "__main__":
    #os.environ["PATH"] = "/usr/local/Caskroom/miniconda/base/envs/GenomeDelta/bin:" + os.environ["PATH"]
    script_path = os.path.abspath("main.sh")

    # Give execute permission and run the script with arguments
    subprocess.call(["chmod", "+x", script_path])

    # Pass command-line arguments to main.sh
    subprocess.call([script_path] + sys.argv[1:])