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
#fortLine="mpif90 -g -Wall -Wno-unused-function \
#-fbounds-check -fcheck=all -mtune=native \
#-fmax-stack-var-size=65536 "  # -O3
fortLine="mpif90 -g -Wall -Wno-unused-function \
-fbounds-check -mtune=native \
-fmax-stack-var-size=65536 "  # -O3
#----------------------------------------------------------------------------------
mkdir -p $binDirLib
mkdir -p $libDirLib
mkdir -p $modDirLib
#----------------------------------------------------------------------------------
# Make mod and binary files
cd $srcDirLib
buildFiles=""
buildFiles=$buildFiles"kinds.f90 "
buildFiles=$buildFiles"rng.f90 "
buildFiles=$buildFiles"logicalMod.f90 "
buildFiles=$buildFiles"keysMod.f90 "
buildFiles=$buildFiles"strings.f90 "
buildFiles=$buildFiles"general.f90 "
buildFiles=$buildFiles"env.f90 "
buildFiles=$buildFiles"arrayFunctions.f90 "
buildFiles=$buildFiles"sort.f90 "
buildFiles=$buildFiles"constants.f90  "
buildFiles=$buildFiles"isotopes.f90 "
buildFiles=$buildFiles"units.f90 "
buildFiles=$buildFiles"mpiSubs.f90 "
buildFiles=$buildFiles"printMod.f90 "
buildFiles=$buildFiles"matrix.f90 "
buildFiles=$buildFiles"basicMaths.f90 "
buildFiles=$buildFiles"rngFunctions.f90 "
buildFiles=$buildFiles"laplaceTransforms.f90 "
buildFiles=$buildFiles"linearAlgebra.f90 "
buildFiles=$buildFiles"calcFunctions.f90  "
buildFiles=$buildFiles"scienceFunctions.f90 "
buildFiles=$buildFiles"solveFunctions.f90 "
buildFiles=$buildFiles"functionPoints.f90 "
buildFiles=$buildFiles"vectors.f90  "
buildFiles=$buildFiles"lmaM.f90 "
buildFiles=$buildFiles"fitting.f90  "
buildFiles=$buildFiles"regression.f90  "
buildFiles=$buildFiles"interpolation.f90 "
buildFiles=$buildFiles"largeInt.f90 "
buildFiles=$buildFiles"newtonGauss.f90 "
buildFiles=$buildFiles"splinesFitting.f90 "
buildFiles=$buildFiles"rngDist.f90 "
buildFiles=$buildFiles"coordFunctions.f90 "
buildFiles=$buildFiles"activityFunctions.f90 "
buildFiles=$buildFiles"plot.f90 "
buildFiles=$buildFiles"geom.f90 "
buildFiles=$buildFiles"potentials.f90 "
buildFiles=$buildFiles"staticCalcs.f90 "
buildFiles=$buildFiles"qe.f90 "
buildFiles=$buildFiles"materials.f90 "
buildFiles=$buildFiles"maths.f90 "
buildFiles=$buildFiles"testMod.f90 "
buildFiles=$buildFiles"libBP.f90 "
commandLine=$fortLine" -J "$modDirLib" -c "
commandLine=$commandLine" "$buildFiles
eval $commandLine
eval "mv "$srcDirLib"/*.o "$binDirLib
# Make library
#commandLine="ar -vr "$libDirLib"/libBP.a "$binDirLib"/*.o "
commandLine="ar -r "$libDirLib"/libBP.a "$binDirLib"/*.o "
echo "Compiling libBP.a:"
echo $commandLine
eval $commandLine
