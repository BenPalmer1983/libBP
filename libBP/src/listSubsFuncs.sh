#!/bin/bash
#################################################################################
# Bash Script
# List subroutines and functions in f90 Modules
#
#################################################################################
# Bash Functions
lc(){
    case "$1" in
        [A-Z])
        n=$(printf "%d" "'$1")
        n=$((n+32))
        printf \\$(printf "%o" "$n")
    esac
}

# Change to working directory (if not already in)
#cd ${0%/*}
# Set useful variables
scriptName="List Fortran Subroutines and Functions"
workingDir=$(pwd)
currentUser=$(whoami)
printHeading=1
# Print
if [ $printHeading == 1 ]; then
  echo "Script:       "$scriptName
  echo "Working dir:  "$workingDir
  echo "User:         "$currentUser
fi
#################################################################################
# Loop through objects in directory
for filename in $workingDir/*.f90; do
  filenameonly=$(basename "$filename")
  filenameshort=${filenameonly:0:25}

  while IFS='' read -r line || [[ -n "$line" ]]; do
    line="$(echo -e "${line}" | sed -e 's/^[[:space:]]*//')"
    line="$(echo -e "${line}" | sed -e 's/[[:space:]]*$//')"
    lineMod=${line:0:6}
    lineModP=${line:0:16}
    lineFunc=${line:0:8}
    lineSub=${line:0:10}
    lineInt=${line:0:9}
    #echo $lineMod
    if [[ $lineMod == "Module" ]]; then
      if [[ $lineModP != "Module Procedure" ]]; then
        echo "========================================================================"
        echo $line"        "$filenameshort
        echo "========================================================================"
      fi
    fi
    if [[ $lineFunc == "Function" ]]; then
      echo "  "$line
    fi
    if [[ $lineSub == "Subroutine" ]]; then
      echo "  "$line
    fi
    if [[ $lineInt == "Interface" ]]; then
      echo "  "$line
    fi


  done < $filename


#  echo $filenameshort" "$ecutwfc" "$ecutrho" "$energy" "$force" "$kpoints" "$smearing"     "$computeTime
done





#################################################################################
