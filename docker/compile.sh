#!/bin/bash

source ~/.bashrc
source /opt/openfoam5/etc/bashrc

cd $HOME/EOF-Library/libs/commSplit
cp 5.0-dev/* .
wmake

cd ../coupleElmer
wmake

if [ "$#" -eq 0 ]; then
    echo "Solver name not provided, using the default one.."
    export SOLVERS="mhdInterFoam"
else
    export SOLVERS=$@
fi

for i in $SOLVERS
do
    cd ../"$i"
    wmake
done
