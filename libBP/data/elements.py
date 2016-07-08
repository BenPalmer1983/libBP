#!/usr/bin/python
#################################################################################
# Python Script
#
#################################################################################
# Import
import os, sys, time, shutil, string;

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

def spacesOnly(inputStr):
# Swaps carriage returns and tabs with spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) not in (9,10,13,32):
      output = output + testChar
    else:
      output = output + " "
  return output

def oneSpace(inputStr):
# Removes spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if i == 1:
      output = output + testChar
    else:
      testCharBefore = subString(inputStr,i-1,1)
      if ord(testChar) == 32:
        if ord(testCharBefore) != 32:
          output = output + testChar
      else:
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



elementsFile = open('elements.txt', 'r')

i=0
fields=[None]*7

for line in elementsFile:
  line = spacesOnly(line)
  line = line.strip()
  line = oneSpace(line)
  if line:
    i = i + 1
    fields = line.split(" ")
    print "      elementsList%atomicNumber("+str(i)+") = "+removeBlanks(fields[0])
    print "      elementsList%atomicSymbol("+str(i)+") = \""+removeBlanks(fields[1])+"\""
    print "      elementsList%elementName("+str(i)+") = \""+removeBlanks(fields[2])+"\""
    print "      elementsList%group("+str(i)+") = "+removeBlanks(fields[3])
    print "      elementsList%period("+str(i)+") = "+removeBlanks(fields[4])
    print "      elementsList%standardAtomicMass("+str(i)+") = "+removeBlanks(fields[5])
    print "      elementsList%density("+str(i)+") = "+removeBlanks(fields[6])
print "      elementsList%elementCount = "+str(i)





#################################################################################
