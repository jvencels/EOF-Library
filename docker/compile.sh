#!/bin/bash

source /opt/openfoam5/etc/bashrc

cd $HOME/EOF-Library/libs/commSplit
cp 5.0-dev/* .
wmake

cd ../coupleElmer
wmake

cd ../mhdInterFoam
wmake
