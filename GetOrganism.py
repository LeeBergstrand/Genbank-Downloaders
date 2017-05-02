#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Description: A simple program that takes a list of genbank geneome accession numbers and finds the name
#			of the taxa that these accession numbers represent.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires the SeqExtract module (included in the Bio-Scripts repository)
#               - All accessions must link to regular nucleotide genome records however if a genome
#                 is shotgun sequenced you must provide the accession to its Whole Genome Shotgun 
#                 Sequence Project record.
#               - Before using this script to access the NCBI's online resources please read the NCBI's 
#                 Entrez User Requirements. If the NCBI finds you are abusing their systems, they can 
#                 and will ban your access! Use the optional email parameter so the NCBI can contact 
#                 you if there is a problem.
#  
# Usage: GetOrganism.py <sequences.txt> [email@mail.com]
# Example: GetOrganism.py mySeqs.txt JBro@YOLO.com
# ----------------------------------------------------------------------------------------
# ===========================================================================================================
# Imports:

import sys

from Bio import Entrez

from SeqExtract import entrezEmail


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
print("Opening sequence list...")
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

for seq in seqList:
    GenbankAccession = seq

    handle = Entrez.esearch(db="nuccore", term=GenbankAccession)
    GenbankSearchResults = Entrez.read(handle)
    if not GenbankSearchResults["IdList"]:
        print(seq + ": Taxa's name not found.")
        continue
    else:
        AccessionGenId = GenbankSearchResults["IdList"][0]

        GenbankInfo = Entrez.esummary(db="nuccore", id=AccessionGenId)
        GenbankSummery = Entrez.read(GenbankInfo)
        TaxaID = GenbankSummery[0]["TaxId"]

        TaxonomyInfo = Entrez.esummary(db="taxonomy", id=TaxaID)
        TaxonomySummary = Entrez.read(TaxonomyInfo)
        TaxaName = TaxonomySummary[0]["ScientificName"]

        print(seq + ": " + TaxaName)
