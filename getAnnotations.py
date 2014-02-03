#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A simple program that takes a list of genbank accesion numbers and downloads 
#           their associated fasta formatted sequences. Stores the sequences in a mutiple
#           sequence fasta file.
#
#           This script requires the Biopython module: http://biopython.org/wiki/Download
#  
# Usage: getGenbankSeqs.py <sequences.txt> <outputName.fasta> <email@mail.com>
# Example: getGenbankSeqs.py mySeqs.txt mySeqs.fasta JBro@YOLO.com
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:
	
import sys
from SeqExtract import entrezEmail
from SeqExtract import getSeqRecords
from SeqExtract import isSSProject
from SeqExtract import extractContigs
from SeqExtract import getProtienAnnotationFasta
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Sequence Downloader"
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + " <sequences.txt> <outputName.fasta> <email@mail.com>\n"
		print "Examples:" + sys.argv[0] + " mySeqs.txt mySeqs.fasta JBro@YOLO.com"
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
#fasta = getProtienAnnotationFasta(seqRecords[0]) # Builds list fasta files.		
for sequence in seqRecords:
	outFile = sequence.id + ".faa"
	try:
		# Attempted to crearte to output file.
		writeFile = open(outFile, "w") 	
		print "Writing " + outFile + " to file..."
		
		#Checks if the accession leads to a WGSS project. 
		#If accession is a WGSS project... 
		if isSSProject(sequence) == True:
			contigList = extractContigs(sequence.id) # Extract all contig accessions. 
			contigRecords = getSeqRecords(contigList) # Extract sequence record object for each contig
			for contig in contigRecords:
				fasta = getProtienAnnotationFasta(contig) # Builds list fasta files.
				for annotation in fasta:
					writeFile.write(annotation)	
		#If accession is a regular genome... 
		else:	
			fasta = getProtienAnnotationFasta(sequence) # Builds list fasta files.
			for annotation in fasta:
				writeFile.write(annotation)		
		writeFile.close()
	except IOError:
		print "Failed to create " + outFile
		exit(1)	
	
print "Done!"





# Attempted to crearte to output file.
#try:
#	writeFile = open(outFile, "w") 	
#	print "Writing sequences to file..."
#	for y in  range(0, len(fasta)):
#		writeFile.write(fasta[y])
#	writeFile.close()
#except IOError:
#	print "Failed to create " + outFile
#	exit(1)	
	
#print "Done!"




	