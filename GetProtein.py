#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Description: A simple program that takes a list of genbank protein accession numbers and downloads
#           their associated sequences. It then stores the genes or complete genomes as a single 
#           peptide fasta file. 
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires the SeqExtract module (included in the Bio-Scripts repository)
#               - All accessions must link to regular nucleotide genbank records (gene or genome),
#                 however if a genome is shotgun sequenced you must provide the accession to its Whole
#                 Genome Shotgun Sequence Project record.
#               - Before using this script to access the NCBI's online resources please read the NCBI's 
#                 Entrez User Requirements. If the NCBI finds you are abusing their systems, they can 
#                 and will ban your access! Use the optional email parameter so the NCBI can contact 
#                 you if there is a problem.
#  
# Usage: GetProtein.py <sequences.faa> [email@mail.com]
# Example: GetProtein.py mySeqs.txt JBro@YOLO.com
# ----------------------------------------------------------------------------------------
# ===========================================================================================================
# Imports:

import sys

from SeqExtract import entrezEmail
from SeqExtract import getSeqRecords


# ===========================================================================================================
# Functions:


# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print("Sequence Downloader")
        print("By Lee Bergstrand\n")
        print("Usage: " + sys.argv[0] + " <sequences.txt> [email@mail.com]")
        print("Examples: " + sys.argv[0] + " mySeqs.txt JBro@YOLO.com\n")
        print("Please Note:")
        print("Before using this script to access the NCBI's online resources please read the NCBI's")
        print("Entrez User Requirements. If the NCBI finds you are abusing their systems, they can")
        print("and will ban your access! Use the optional email parameter so the NCBI can contact")
        print("you if there is a problem.")
        sys.exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# ===========================================================================================================
# Main program code:

# House keeping...
argsCheck(3)  # Checks if the number of arguments are correct.
entrezEmail(sys.argv[2])  # Sets up arguments email require for genbank file extraction.

# Stores file one for input checking.
print(">> Opening sequence list...")
inFile = sys.argv[1]

# File extension check
if not inFile.endswith(".txt"):
    print("[Warning] " + inFile + " may not be a txt file!")

# Reads sequence file list and stores it as a string object. Safely closes file:
try:
    with open(inFile, "r") as newFile:
        sequences = newFile.read()
        newFile.close()
except IOError:
    print("Failed to open " + inFile)
    sys.exit(1)

seqList = sequences.splitlines()  # Splits string into a list. Each element is a single line from the string.

print("You have listed", len(seqList), "sequences. They are:")
print(sequences + "\n\n")

seqRecords = getSeqRecords(seqList, database_type="protein")  # Gets sequence record objects from NCBI using the sequence list as reference.

outFile = inFile + ".faa"
try:
    # Attempted to create to output file.
    writeFile = open(outFile, "w")
    print("Writing " + outFile + " to file...")
    for sequence in seqRecords:
        writeFile.write(sequence.format("fasta"))  # Write genome as fasta to file.
    writeFile.close()
except IOError:
    print("Failed to create " + outFile)
    sys.exit(1)

print("Done!")
