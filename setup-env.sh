#!/bin/bash

ENV_NAME="xmm_epic_env"
YAML_FILE="environment.yml"

# Try Python 3.13 first
cat > $YAML_FILE <<EOF
name: $ENV_NAME
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.13
  - numpy
  - pandas
  - matplotlib
  - astropy
  - requests
  - wget
  - pip
  - pip:
      - glob2
EOF

echo "Trying to create environment with Python 3.13..."
conda env create -f $YAML_FILE

if [ $? -eq 0 ]; then
    echo "Environment created with Python 3.13."
    conda activate $ENV_NAME
    # Try importing heasarc/xmmsas or check their install
    # python -c "import heasarc"  # Uncomment if you have a test import
    # If import fails, set status=1
    status=0
else
    status=1
fi

if [ $status -ne 0 ]; then
    echo "Python 3.13 or dependencies not compatible. Falling back to Python 3.9..."
    cat > $YAML_FILE <<EOF
name: $ENV_NAME
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.9
  - numpy
  - pandas
  - matplotlib
  - astropy
  - requests
  - wget
  - pip
  - pip:
      - glob2
EOF
    conda env create -f $YAML_FILE
    if [ $? -eq 0 ]; then
        echo "Environment created with Python 3.9."
    else
        echo "Failed to create environment with both Python 3.13 and 3.9."
        exit 1
    fi
fi

echo "Done. Activate with: conda activate $ENV_NAME"