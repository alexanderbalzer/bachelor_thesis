#!/bin/bash

# Set the environment YAML file and the Python script
ENV_FILE="bachelor_env.yml"
PYTHON_SCRIPT="main.py"
TAR_FILE="MitoFates_1.2(1).tar.gz"
EXTRACT_DIR="MitoFates"

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

# Unpack MitoFates
if [ ! -d "$EXTRACT_DIR" ]; then
    echo "Extracting $TAR_FILE..."
    mkdir -p "$EXTRACT_DIR"
    tar -xzf "$TAR_FILE" -C "$EXTRACT_DIR" --strip-components=1
else
    echo "$EXTRACT_DIR already exists. Skipping extraction."
fi

# Run the script
python "$PYTHON_SCRIPT"


# Determine directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Navigate to target folder: parent of script dir, then into a child folder (e.g., "analysis")
TARGET_DIR="$SCRIPT_DIR/../machine_learning"
cd "$TARGET_DIR" || { echo "Target directory not found: $TARGET_DIR"; exit 1; }

# Run four Python scripts
python create_feature_matrix.py
python supervised_machine_learning_per_GO_term.py
python logreg_results_to_barchart.py
python logreg_results_to_heatmap.py
