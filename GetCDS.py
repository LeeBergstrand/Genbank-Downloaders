#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Description: A simple program that takes a list of nucleotide genbank accession numbers and
#			downloads the Coding Sequences (CDS) contained within the sequences linked by  
#  			those accessions. It then stores these CDSs within a multi-sequence protein fasta
#			file. The script also creates a CSV file containing some essential info about each CDS.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires the SeqExtract module (included in the Genbank-Downloaders repository)
#               - All accessions must link to regular nucleotide genbank records (gene or genome),
#                 however if a genome is shotgun sequenced you must provide the accession to its Whole
#                 Genome Shotgun Sequence Project record.
#               - Before using this script to access the NCBI's online resources please read the NCBI's 
#                 Entrez User Requirements. If the NCBI finds you are abusing their systems, they can 
#                 and will ban your access! Use the optional email parameter so the NCBI can contact 
#                 you if there is a problem.
#  
# Usage: Get16S.py <sequences.txt> [email@mail.com]
# Example:Get16S.py mySeqs.txt JBro@YOLO.com
# ----------------------------------------------------------------------------------------
# ===========================================================================================================
# Imports:
	
import sys
import csv
from SeqExtract import entrezEmail
from SeqExtract import getSeqRecords
from SeqExtract import isSSProject
from SeqExtract import extractContigs
from SeqExtract import getProtienAnnotationFasta
from SeqExtract import getProtienAnnotationCSV
# ===========================================================================================================
# Functions:


# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Coding Sequence Downloader"
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + " <sequences.txt> [email@mail.com]"
		print "Examples: " + sys.argv[0] + " mySeqs.txt JBro@YOLO.com\n"
		print "Please Note:"
		print "Before using this script to access the NCBI's online resources please read the NCBI's" 
		print "Entrez User Requirements. If the NCBI finds you are abusing their systems, they can" 
		print "and will ban your access! Use the optional email parameter so the NCBI can contact" 
		print "you if there is a problem."
		sys.exit(1) # Aborts program. (exit(1) indicates that an error occurred)
# ===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(3) # Checks if the number of arguments are correct.
entrezEmail(sys.argv[2]) # Sets up arguments email require for genbank file extraction.
	
# Stores file one for input checking.
print ">> Opening sequence list..."
inFile = sys.argv[1]

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
	sys.exit(1)

seqList = sequences.splitlines()  # Splits string into a list. Each element is a single line from the string.

print "You have listed", len(seqList), "sequences. They are:"
print sequences + "\n\n"
	
seqRecords = getSeqRecords(seqList)  # Gets sequence record objects from NCBI using the sequence list as reference.
		
for sequence in seqRecords:
	outFile = sequence.id + ".faa"
	outCSV = sequence.id + ".csv"
	try:
		# Attempted to create to output file.
		writeFile = open(outFile, "w") 	
		print "Writing " + outFile + " to file..."
		csvFile = open(outCSV, "w") 	
		CSVWriter = csv.writer(csvFile)
		print "Writing " + outCSV + " to file..."
		
		# Checks if the accession leads to a WGSS project.
		# If accession is a WGSS project...
		if isSSProject(sequence):
			contigList = extractContigs(sequence.id)  # Extract all contig accessions.
			contigRecords = getSeqRecords(contigList)  # Extract sequence record object for each contig.
			for contig in contigRecords:
				fasta = getProtienAnnotationFasta(contig)  # Builds list fasta files.
				csvRows = getProtienAnnotationCSV(contig)  # Builds list of csv rows.
				for annotation in fasta:
					writeFile.write(annotation)	
				for row in csvRows:
					CSVWriter.writerow(row)
		# If accession is a regular genome...
		else:
			OrganismGenomeLength = len(sequence.seq)  # Gets Genome Length
			fasta = getProtienAnnotationFasta(sequence)  # Builds list fasta files.
			csvRows = getProtienAnnotationCSV(sequence)  # Builds list of csv rows.
			for annotation in fasta:
				writeFile.write(annotation)	
			for row in csvRows:
				row.append(OrganismGenomeLength)  # Appends genome length to the rest of the csv file.
				CSVWriter.writerow(row)	
		writeFile.close()
		csvFile.close()
	except IOError:
		print "Failed to create " + outFile
		sys.exit(1)	
print "Done!"