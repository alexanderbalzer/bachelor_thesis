#!/bin/bash

# Set the environment YAML file and the Python script
ENV_FILE="bachelor_env.yml"
PYTHON_SCRIPT="main.py"

# Extract environment name from YAML (assumes first line is `name: env_name`)
ENV_NAME=$(grep '^name:' "$ENV_FILE" | cut -d' ' -f2)

# Check if environment already exists
if ! conda info --envs | grep -q "^$ENV_NAME[[:space:]]"; then
    echo "Creating Conda environment '$ENV_NAME' from $ENV_FILE..."
    conda env create -f "$ENV_FILE"
else
    echo "Conda environment '$ENV_NAME' already exists."
fi

# Activate the environment and run the Python script
echo "Activating environment '$ENV_NAME' and running '$PYTHON_SCRIPT'..."
# Make sure conda is available in this shell
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Run the script
python "$PYTHON_SCRIPT"
