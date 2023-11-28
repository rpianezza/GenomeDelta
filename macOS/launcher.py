import os
import subprocess
import sys

if __name__ == "__main__":

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to "main.sh"
    script_path = os.path.join(script_dir, "main.sh")

    # Give execute permission and run the script with arguments
    subprocess.call(["chmod", "+x", script_path])

    # Pass command-line arguments to main.sh
    subprocess.call([script_path] + sys.argv[1:])