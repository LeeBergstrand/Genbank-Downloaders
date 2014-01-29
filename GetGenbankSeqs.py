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
#==============================================================================================================================
#Imports:
	
import sys
from Bio import SeqIO
from Bio import Entrez
#==============================================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Sequence Downloader"
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + " <sequences.txt> <outputName.fasta> <email@mail.com>\n"
		print "Examples:" + sys.argv[0] + " mySeqs.txt mySeqs.fasta JBro@YOLO.com"
		exit(1) # Aborts program. (exit(1) indicates that an error occured)

# 2: When passed an array of accessions of NCBI returns a list of sequence objects matching those accessions.
def getSeqRecords(seqList):
	try: 
		print "Requesting sequence data from genbank..."
		handle = Entrez.efetch(db="protein", id=seqList, rettype="gb", retmode="genbank") # Gets genbank files and stores them.
		print "Starting download..."
		SeqRecords=list(SeqIO.parse(handle,"genbank")) # Creates a list of SeqRecord objects from genbank files stored in handle.
		print "Download Complete."
		handle.close() # Closes handle since it is no longer needed.
	except IOError:
		print "Failed to connect to NCBI server. "
		exit(1)
		
	return SeqRecords
#==============================================================================================================================
# Main program code:
	
# House keeping:
argsCheck(4)
	
# Stores file one for input checking.
print ">> Opening sequence list..."
inFile  = sys.argv[1]
outFile = sys.argv[2]
Entrez.email = sys.argv[3]

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
	
try:
	writeFile = open(outFile, "w") 	
except IOError:
		print "Failed to create " + outFile
		exit(1)
		
SeqRecords = getSeqRecords(seqList)
		
print "Writing sequences to file..."
try:
	for x in range(0, len(SeqRecords)): # For each sequence
		features = SeqRecords[x].features # Each sequence has a list (called features) that stores seqFeature objects.
		for y in range(0, len(features)): # For each feature on the sequence
			if features[y].type == "CDS": # CDS means coding sequence (These are the only feature we're interested in)
				featQualifers = features[y].qualifiers # Each feature contains a dictionary called quailifiers which contains data about         
				                                       # the sequence feature (for example the translation)
				
				# Gets the required qualifers. Uses featQualifers.get to return the quatifer or a default value if the quatifer				# is not found. Calls strip to remove unwanted brackets and ' from quantifer before storing it as a string.
				protein_id = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')
				if protein_id == 'no_protein_id':
					continue # Skips the iteration if protien has no id.
				gene = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
				product = str(featQualifers.get('product','no_product_name')).strip('\'[]')
				translated_protein = str(featQualifers.get('translation','no_translation')).strip('\'[]')
				
				fasta = ">" + protein_id + " " + gene + "-" + product + "\n" + translated_protein + "\n"
				writeFile.write(fasta)
	writeFile.close()
	print "File writing complete."
except IOError:
	print "Failed to write to output file " + outFile
	exit(1)	
	
print "Done!"


	