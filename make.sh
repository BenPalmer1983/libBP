#!/bin/bash
cd ${0%/*}
workingDir=$(pwd)
#src
srcDir=$workingDir"/src"
#Directories
modDir=$workingDir"/mod"
binDir=$workingDir"/bin"
mkdir -p $modDir
mkdir -p $binDir
#-----------Library Start--------
libDir=$workingDir"/libBP"
modDirLib=$libDir"/mod"
libDirLib=$libDir"/lib"
libLink="-L"$libDirLib" -lBP"
#Include library
"$workingDir/makeLib.sh"
#-----------Library End--------
fortLine="mpif90 -O3 -g -Wall -Wno-unused-function \
-fbounds-check -fcheck=all -mtune=native "
cd $srcDir
commandLine=$fortLine" -J "$modDir"  -I"$modDirLib" testProg.f90 "$libLink" -o "$binDir"/test.x"
eval $commandLine

