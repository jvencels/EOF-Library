#!/bin/bash

source /opt/openfoam5/etc/bashrc
source ~/.bashrc

containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

EOFdir='..' # EOF-Library folder

while getopts "f:v:s:" flag; do
  case "$flag" in
    f  ) EOFdir=$OPTARG ;;
    v  ) OFvers=$OPTARG ;;
    s  ) SOLVERS+=("$OPTARG") ;;
    \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
    :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    *  ) echo "Unexpected option ${flag}. Valid options are -f, -v, -s" && exit 1;;
  esac
done

validOFvers=("2.4.0" "3.0.1" "4.1" "5.0-stable" "5.0-dev")

cd $EOFdir/libs/commSplit
echo | pwd

if containsElement "$OFvers" "${validOFvers[@]}"; then
  echo "Configuring EOF-Library for OpenFOAM version: $OFvers"
  cp $OFvers/* .
  wclean && wmake
elif [ "$OFvers" == "dev" ]; then
  echo "Configuring EOF-Library for OpenFOAM version: dev"
else
	echo "ERROR: OpenFOAM version $OFvers not supported. Valid versions are:"
    printf "%s\n" "${validOFvers[@]}"
    printf "dev\n"
    exit 1
fi

echo "Compiling coupler.."
cd ../coupleElmer
wclean && wmake

echo "Compiling solvers.."
for i in ${SOLVERS[@]}
do
  if test -d ../$i; then
    echo "$i"
    cd ../$i
    wclean && wmake
  else
    echo "ERROR: Solver $i does not exist!"
    exit 1
  fi
done
