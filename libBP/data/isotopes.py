#!/usr/bin/python
#################################################################################
# Python Script
#
#################################################################################
# Import
import os, sys, time, shutil;

def subString(inputStr, start, length):
  a = start-1
  output = inputStr[a:a+length]
  return output

def removeSpaces(inputStr):
# Removes spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) != 32:
      output = output + testChar
  return output

def removeBlanks(inputStr):
# Removes spaces, tabs, returns
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) not in (9,10,13,32):
      output = output + testChar
  return output

def cleanReal(inputStr):
# Removes spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) == 43:
      output = output + testChar
    if ord(testChar) == 45:
      output = output + testChar
    if ord(testChar) == 46:
      output = output + testChar
    if ord(testChar) >= 48 and ord(testChar) <= 57 :
      output = output + testChar
  return output



isotopeFile = open('isotopes.txt', 'r')

i=0
atomicNumber=[None]*50000
atomicSymbol=[None]*50000
massNumber=[None]*50000
relativeAtomicMass=[None]*50000
isotopicComposition=[None]*50000
standardAtomicWeight=[None]*50000



for line in isotopeFile:
  if(line[0:15]=="Atomic Number ="):
    i = i + 1
    atomicNumber[i] = cleanReal(removeSpaces(line[16:20]))
    atomicNumber[i] = atomicNumber[i].upper()
  if(line[0:15]=="Atomic Symbol ="):
    tempSymbol = removeBlanks(line[16:20])
    atomicSymbol[i] = tempSymbol[0:12]
  if(line[0:13]=="Mass Number ="):
    massNumber[i] = removeSpaces(line[13:17])
  if(line[0:22]=="Relative Atomic Mass ="):
    relativeAtomicMass[i] = cleanReal(line[23:40])
    if(relativeAtomicMass[i]==""):
      relativeAtomicMass[i] = massNumber[i]
  if(line[0:22]=="Isotopic Composition ="):
    isotopicComposition[i] = cleanReal(line[23:40])
    if(isotopicComposition[i]==""):
      isotopicComposition[i] = 0.0
  if(line[0:24]=="Standard Atomic Weight ="):
    standardAtomicWeight[i] = cleanReal(line[25:40])


print ""

lineCountMax = i
lineCount = str(i)
print "      isotopesList%isotopeCount = "+lineCount

i=0
while i<lineCountMax:
  i = i + 1
  print "      isotopesList%atomicSymbol("+str(i)+") = \""+removeBlanks(atomicSymbol[i])+"\""
  print "      isotopesList%atomicNumber("+str(i)+") = "+removeBlanks(atomicNumber[i])
  print "      isotopesList%massNumber("+str(i)+") = "+removeBlanks(massNumber[i])
  print "      isotopesList%isotopicComposition("+str(i)+") = "+removeBlanks(str(isotopicComposition[i]))
  print "      isotopesList%relativeAtomicMass("+str(i)+") = "+removeBlanks(str(relativeAtomicMass[i]))
  print "      isotopesList%standardAtomicMass("+str(i)+") = "+removeBlanks(str(standardAtomicWeight[i]))





#################################################################################
