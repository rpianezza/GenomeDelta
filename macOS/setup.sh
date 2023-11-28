#!/bin/bash

# Download the GitHub repository
git clone https://github.com/rpianezza/GenomeDelta.git

# Move to the downloaded folder
cd GenomeDelta/macOS

# Create the conda environment
conda env create -f set-env.yml

# Activate the conda environment
conda activate GenomeDelta

# Get the path to the conda environment's bin directory
conda_bin_path=$(conda info --envs | grep "^GenomeDelta " | awk '{print $NF}')/bin

# Get the path to "main.sh"
main_script_path=$(pwd)/main.sh

# Create the executable file
python -m pip install pyinstaller
pyinstaller --onefile --name GenomeDelta launcher.py -- --main "$main_script_path"

# Copy the executable file to the conda environment's bin directory
mv dist/GenomeDelta "$conda_bin_path"

# Print information about the setup
echo "Setup complete. Type "conda activate GenomeDelta" to activate the environment."
echo "To launch the software, type GenomeDelta in the environment and add the requested arguments (--fq/--bam, --fa, --of, --t)"