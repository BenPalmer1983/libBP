#!/bin/bash
# ----------------------------------
# University of Birmingham
# Ben Palmer
# ----------------------------------
cd ${0%/*}
workingDirLib=$(pwd)
srcDirLib=$workingDirLib"/libBP/src"
modDirLib=$workingDirLib"/libBP/mod"
libDirLib=$workingDirLib"/libBP/lib"
binDirLib=$workingDirLib"/libBP/bin"
fortLine="mpif90 -O3 -g -Wall -Wno-unused-function \
-fbounds-check -fcheck=all -mtune=native "  # -O3
#----------------------------------------------------------------------------------
mkdir -p $binDirLib
mkdir -p $libDirLib
mkdir -p $modDirLib
#----------------------------------------------------------------------------------
# Make mod and binary files
cd $srcDirLib
commandLine=$fortLine" -J "$modDirLib" -c "
commandLine=$commandLine"kinds.f90 strings.f90 general.f90 "
commandLine=$commandLine"arrayFunctions.f90 constants.f90 "
commandLine=$commandLine"units.f90 printMod.f90 matrix.f90 basicMaths.f90 rng.f90 "
commandLine=$commandLine"laplaceTransforms.f90 linearAlgebra.f90 calcFunctions.f90 "
commandLine=$commandLine"solveFunctions.f90 functionPoints.f90 vectors.f90 "
commandLine=$commandLine"lmaM.f90 regression.f90 "
commandLine=$commandLine"interpolation.f90 newtonGauss.f90 "
commandLine=$commandLine"fitting.f90 rngDist.f90 coordFunctions.f90 "
commandLine=$commandLine"activityFunctions.f90 "
commandLine=$commandLine"splines.f90 fitting.f90 plot.f90 geom.f90 "
commandLine=$commandLine"maths.f90 libBP.f90 "
eval $commandLine
eval "mv "$srcDirLib"/*.o "$binDirLib
# Make library
commandLine="ar -vr "$libDirLib"/libBP.a "$binDirLib"/*.o " 
eval $commandLine 