#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A simple program that takes a list of genbank genome accesion numbers and downloads 
#           the geneome as fasta formatted sequence. For shotgun sequenced data  it stores the 
#           contigs in a mutiple sequence fasta file.
#
#           This script requires the Biopython module: http://biopython.org/wiki/Download
#  
# Usage: getGenbankSeqs.py <sequences.txt> <email@mail.com>
# Example: getGenbankSeqs.py mySeqs.txt JBro@YOLO.com
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:
	
import sys
from SeqExtract import entrezEmail
from SeqExtract import getSeqRecords
from SeqExtract import getProtienAnnotationFasta
from SeqExtract import isSSProject
from SeqExtract import extractContigs
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Sequence Downloader"
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + " <sequences.txt> <email@mail.com>\n"
		print "Examples:" + sys.argv[0] + " mySeqs.txt JBro@YOLO.com"
		exit(1) # Aborts program. (exit(1) indicates that an error occured)
#===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(3) # Checks if the number of arguments are correct.
entrezEmail(sys.argv[2]) # Sets up arguments email require for genbank file extraction.
	
# Stores file one for input checking.
print ">> Opening sequence list..."
inFile  = sys.argv[1]

# File extension check
if not inFile.endswith(".txt"):
	print "[Warning] " + inFile + " may not be a txt file!"

# Reads sequence file list and stores it as a string object. Safely closes file.try:
try:	
	with open(inFile,"r") as newFile:
		sequences = newFile.read()
		newFile.close()
except IOError:
	print "Failed to open " + inFile
	exit(1)

seqList = sequences.splitlines() # Splits string into a list. Each element is a single line from the string.

print "You have listed", len(seqList), "sequences. They are:"
print sequences + "\n\n"
	
seqRecords = getSeqRecords(seqList) # Aquires list of sequence record objects from NCBI us sequence list as reference.

for sequence in seqRecords:
	outFile = sequence.id + ".fna"
	try:
		# Attempted to crearte to output file.
		writeFile = open(outFile, "w") 	
		print "Writing " + outFile + " to file..."
		
		# Checks if the accession leads to a WGSS project. 
		# If accession is a WGSS project... 
		if isSSProject(sequence) == True:
			contigList = extractContigs(sequence.id) # Extract all contig accessions. 
			contigRecords = getSeqRecords(contigList) # Extract sequence record object for each contig
			for contig in contigRecords:
				writeFile.write(contig.format("fasta")) # Write each contig to the same file in fasta format.
		# If accession is a regular genome... 
		else:	
			writeFile.write(sequence.format("fasta")) # Write genome as fasta to file
		writeFile.close()
	except IOError:
		print "Failed to create " + outFile
		exit(1)	
	
print "Done!"




	