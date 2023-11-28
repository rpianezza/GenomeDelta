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

# Create the executable file
python -m pip install pyinstaller
pyinstaller launcher.py --onefile --name GenomeDelta

# Copy the executable file to the conda environment's bin directory
mv dist/GenomeDelta "$conda_bin_path"

# Print information about the setup
echo "Setup complete."
echo "To launch the software, type GenomeDelta and add the requested arguments (--fq/--bam, --fa, --of, --t)"