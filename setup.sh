#!/bin/bash

mkdir -p ~/eof_runs

# Save the current directory as EOF_HOME
EOF_HOME=$(pwd)

# Add environment variables to .bashrc
echo "export EOF_RUNS_DIR=~/eof_runs" >> ~/.bashrc
echo "export EOF_HOME=$EOF_HOME" >> ~/.bashrc

# Alias to build the local EOF Docker image
echo "alias eof_build='cd \$EOF_HOME && docker build --network host --platform linux/amd64 -f docker/Dockerfile.eof_local -t eof_local .'" >> ~/.bashrc

# Alias to run the EOF container interactively
echo "alias eof_run='docker run --rm -it --platform linux/amd64 -e HOST_USER_ID=\$(id -u) -e HOST_USER_GID=\$(id -g) -v \$EOF_RUNS_DIR:/home/openfoam/EOF-Library/runs eof_local'" >> ~/.bashrc

# The master command to build and run the EOF container
echo "alias eof_start='eof_build && eof_run'" >> ~/.bashrc

echo "Setup complete! Please run: source ~/.bashrc"