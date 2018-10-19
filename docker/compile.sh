#!/bin/bash

containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

# Run from main directory
cd "$(dirname "$0")"/..

while getopts ":v:s:" flag; do
  case "$flag" in
    v  ) OFvers=$OPTARG ;;
    s  ) SOLVERS+=("$OPTARG") ;;
    \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
    :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    *  ) echo "Unexpected option ${flag}. Valid options are -v, -s" && exit 1;;
  esac
done

validOFvers=("2.4.0" "3.0.1" "4.1" "5.0" "5.0-dev" "6" "dev" "")

echo | pwd

if containsElement "$OFvers" "${validOFvers[@]}"; then
  echo "Configuring EOF-Library for OpenFOAM version: $OFvers"
  wclean
  if [ "$OFvers" = "6" ] || [ "$OFvers" = "dev" ]; then
    continue
  elif [ "$OFvers" = "5.0-dev" ]; then
    echo "Applying OF5DEV patch!"
    export OF5DEV=OF5DEV
  fi
  wmake
else
  echo "ERROR: OpenFOAM version $OFvers not supported. Valid versions are:"
  printf "%s\n" "${validOFvers[@]}"
  exit 1
fi

echo "Compiling coupler.."
cd ../coupleElmer
wclean && wmake

echo "Compiling solvers.."
cd ../../solvers
for i in ${SOLVERS[@]}
do
  if test -d $i; then
    echo "$i"
    wclean $i && wmake $i
  else
    echo "ERROR: Solver $i does not exist!"
    exit 1
  fi
done
