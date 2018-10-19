#!/bin/bash

# Script for modifying Dockerfiles and building them

# OPTIONS:
# -b  = build type "local" or "pull"
# -f  = path to Dockerfile
# -t  = name and tag for Docker image

# FOR TESTING:
# ./dockerBuild.sh -b local -f Dockerfile_eof_of6_debug -t eof_of6:debug

# FOR DEPLOYING:
# ./dockerBuild.sh -b pull -f Dockerfile_eof_of6 -t eof_of6:latest
# ./dockerBuild.sh -b pull -f Dockerfile_eof_swak4foam_of6 -t eof_swak4foam:of6

# Run from main directory
cd "$(dirname "$0")"/..

# Function for searching string in string array
containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

# Define options
while getopts ":b:f:t:" flag; do
  case "$flag" in
    b  ) buildType=$OPTARG ;;  # Build type
    f  ) dockerFile=$OPTARG ;; # DockerFile
    t  ) nameTag=$OPTARG ;;   # Name and tag
    \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
    :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    *  ) echo "Unexpected option ${flag}. Valid options are -f, -v, -s" && exit 1;;
  esac
done

# Check if Dockerfile exists
if [ ! -f $dockerFile ]; then
  if [ ! -f docker/$dockerFile ]; then
    echo "ERROR: File '$dockerFile' not found!"
    exit 1
  else
    dockerFile="docker/$dockerFile"
  fi
fi

# Local EOF-Library or pull from  github
validBuildTypes=("local" "pull")
if containsElement "$buildType" "${validBuildTypes[@]}"; then
  echo "Docker build type: $buildType"
  if [ "$buildType" = "local" ]; then
    sed 's@RUN git clone https://github.com/jvencels/EOF-Library.git@ADD ./ EOF-Library@g' $dockerFile > $dockerFile-tmp
  fi
else
  echo "ERROR: Docker build type '$buildType' not supported. Valid versions are:"
  printf "%s\n" "${validBuildTypes[@]}"
  exit 1
fi

# Check name and tag
if [ -z "$nameTag" ]
then
  echo "ERROR: Docker name:tag is not specified! Tag is optional.."
  exit 1
fi

echo "Building $dockerFile-tmp"
docker build --no-cache -t eof-library/$nameTag -f $dockerFile-tmp .
rm $dockerFile-tmp
