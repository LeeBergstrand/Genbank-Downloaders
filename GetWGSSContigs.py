#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: 
#
#           This script requires the Biopython module: http://biopython.org/wiki/Download
#  
# Usage: 
# Example:
#----------------------------------------------------------------------------------------
import sys
import re
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "" # Temporary, script will request user email in later version.

# Regexs for contig accession extraction
WGSSRangeRegex = re.compile("/wgs=\['{1}[a-zA-Z0-9]{12,14}',[ ]{1}'[a-zA-Z0-9]{12,14}']{1}")
WGSSRangeAccessionsRegex = re.compile("[a-zA-Z0-9]{12,14}")
AccessionBaseRegex = re.compile("^[a-zA-Z]{4}\d{2}")

WGSSProjectId = ["CBMO000000000"]
WGSSProjectId.append("BAEC00000000")

try: 
	print "Requesting sequence data from genbank..."
	handle = Entrez.efetch(db="protein", id = WGSSProjectId, rettype="gb", retmode="genbank") # Gets genbank files and stores them.
	print "Starting download..."
	SeqRecords=list(SeqIO.parse(handle,"genbank")) # Creates a list of SeqRecord objects from genbank files stored in handle.
	print "Download Complete."
	handle.close() # Closes handle since it is no longer needed.
except IOError:
	print "Failed to connect to NCBI server. "
	exit(1)
	
for WGSS in SeqRecords:
	
	# Takes WGSS seqRecord object and converts it to a string. Uses WGSRangeRegex's findAll method to return  
	# a list containing the WGS contig accession range (should be one occurence of this). Pulls out the WGS 
	# contig accession range, as a string, from this list (first index). Applys WGSSRangeAccessionsRegex's 
	# findAll method to return a list with both the max and the min accession for the WGS contig range. 
	WGSSRange = WGSSRangeAccessionsRegex.findall(WGSSRangeRegex.findall(str(WGSS))[0])
	
	AccessionBase = (AccessionBaseRegex.findall(WGSSRange[0])[0]) # Extracts the accession base from first contig accession number.
	
	# Takes both the the min and max accession and slices off (using python's string slice syntax s[start:end:step]) 
	# the accession base code leaving the numerical difference between the contigs. Converts these differences to integers.
	WGSSRangeMin = int(WGSSRange[0][6:]) 
	WGSSRangeMax = int(WGSSRange[1][6:])
	
	# WGSS accession number length actually varies. Its normally 12 characters but I have seen 13 before. 
	# The code block below accounts for this.
	zeroOffset = 6
	accessionLength = len(WGSSRange[0])
	if accessionLength != 12:
		zeroOffset = accessionLength - 6 # 6 is the length of the standard accession base.
	
	contigList = []
	
	for x in range(WGSSRangeMin, (WGSSRangeMax + 1)):
		contigAccession = AccessionBase
		contigAccession += ("{0:0" + str(zeroOffset) + "d}").format(x) # Uses zero offset to make accessions proper length.
		contigList.append(contigAccession)
	
	print contigList

	


	

