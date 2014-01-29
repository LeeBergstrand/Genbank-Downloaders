#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A simple script that modifies the elements inside a column of a CSV by
#           using a regular expression to find and replace charaters in those elements.
#
# Usage: CSVmod.py <input.csv> <output.csv> <columnNumber> <regex> <replace>
# Example: CSVmod.py myInput.csv myOutput.csv 6 ^[\t]+|[\t]$ replacement
#----------------------------------------------------------------------------------------

import csv
import sys
import re

# If in proper number of arguments are passed gives instructions on proper use.
if len(sys.argv) < 6 or len(sys.argv) > 6:
	print "CSV Modifier"
	print "By Lee Bergstrand\n"
	print "A simple script that modifies the elements inside a column of a CSV by\nusing a regular expression to find and replace charaters in those elements.\n"
	print "Usage: " + sys.argv[0] + " <input.csv> <ouput.csv> <columnNumber> <regex> <replacement>"
	print "Examples: " + sys.argv[0] + ' myInput.csv myOutput.csv 6 ^[dec]+|[6]+ replacement'
	exit(1) # Aborts program. (exit(1) indicates that an error occured)

# Stores stores argument data.
print ">> Opening CSV file..."
inFile  = sys.argv[1]
outFile = sys.argv[2]
column  = (int(sys.argv[3]) - 1) # Converts user specified column number (eg. 1-10) to proper list index
regex   = sys.argv[4]
replace = sys.argv[5]

# File extension check
if not inFile.endswith(".csv"):
	print "[Warning] " + inFile + " may not be a CSV file!"
else:
	print ">> Good file extention."	
	
# Opens CSV file for reading.
try:
	readFile = open(inFile, "r")
	reader  = csv.reader(readFile) # opens file with csv module which takes into account verying csv formats and parses correctly
	print ">> Good CSV file."
except IOError:
	print "Failed to open " + inFile
	exit(1)

# Determines if user specifies a column that does not exist and if he/she does exits
if (len(next(reader)) < (column + 1)) or (column < 0): # "column + 1" converts list index back to number of columns 
	print "Error: The column you want to change does not exist."
	print "Exiting..."
	exit(1)
	
readFile.seek(0) # Resets file iterator back to 0 so program does not skip first row.

# Opens CSV file for writing.
try:
	writeFile = open(outFile, "w") 	
	writer = csv.writer(writeFile)
	print ">> Output file created."
	print ">> Writing Data..."
except IOError:
	print "Failed to create " + outFile
	exit(1)

# Modifies changes elements in column according to regex. New modified rows to new file.
for row in reader:
	newElement = re.sub(regex, replace, row[column]) # Uses regex to replace a portion of an element string with a replacement string.
	row[column] = newElement
	writer.writerow(row) # writes new row to outgoing file.

print ">> Finished Writing."
print ">> Closing..."		
readFile.close()
writeFile.close()

print ">> Done!"
		
