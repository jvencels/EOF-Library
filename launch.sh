#!/bin/bash

# Go to the directory that this script is in
cd "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")"

# Build, and exit strictly on failure
make build || exit 1

docker run --rm -it --platform linux/amd64 -e HOST_USER_ID=$(id -u) -e HOST_USER_GID=$(id -g) -v ./out/:/home/openfoam/EOF-Library/runs eof_local:latest
# docker run --rm -it --platform linux/amd64 -e HOST_USER_ID=$(id -u) -e HOST_USER_GID=$(id -g) eof_local:latest

