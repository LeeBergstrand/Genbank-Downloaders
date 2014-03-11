#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A simple program that takes a nucleotide FASTA file and returns the exact same FASTA
#			file with a reverse complemented sequence. Also works with multisequence FASTA files.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#     
# Usage: getGenbankSeqs.py <sequences.fna> 
# Example: getGenbankSeqs.py mySeqs.fna
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:
	
import sys
from Bio import SeqIO
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Takes a nucleotide FASTA file and returns the exact same FASTA file with a reverse complemented sequence."
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + " <sequences.txt> [email@mail.com]"
		print "Examples: " + sys.argv[0] + " mySeqs.txt JBro@YOLO.com\n"
		exit(1) # Aborts program. (exit(1) indicates that an error occurred)
#===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(2) # Checks if the number of arguments are correct.
	
# Stores file one for input checking.
print ">> Opening FASTA file..."
inFile  = sys.argv[1]
outFile = inFile + ".out"

# File extension check
if not inFile.endswith(".fna"):
	print "[Warning] " + inFile + " may not be a FASTA file!"
	
try:
	writer = open(outFile, "w")
	handle = open(inFile, "rU")
	SeqRecords = SeqIO.parse(handle, "fasta")
	print ">> Creating Reverse Complement..."
	for record in SeqRecords:
		reverseCompSeq = record.seq.reverse_complement()
		record.seq = reverseCompSeq
		writer.write(record.format("fasta"))
	handle.close()
	writer.close()
except IOError:
	print "Failed to open " + inFile + " or " + outFile
	exit(1)
	
print ">> Done."