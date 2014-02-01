#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A module that contains functions that extract and modify genbanks records. 
#
#           This module requires the Biopython module: http://biopython.org/wiki/Download
#------------------------------------------------------------------------------------------------------------
#============================================================================================================
# Imports:

import re
from Bio import SeqIO
from Bio import Entrez
#============================================================================================================
# Functions:
	
# 1: Set up email.
def entrezEmail(email):
	Entrez.email = email	
#------------------------------------------------------------------------------------------------------------
				
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
#------------------------------------------------------------------------------------------------------------
	
# 3: When passed an array of sequence record objects returns an array of fasta strings for each annotation.
def getProtienAnnotationFasta(SeqRecords):
	fasta = []
	# Loops through sequence records and extracts required information about protein annotions.
	for x in range(0, len(SeqRecords)): # For each sequence
		features = SeqRecords[x].features # Each sequence has a list (called features) that stores seqFeature objects.
		for y in range(0, len(features)): # For each feature on the sequence
			if features[y].type == "CDS": # CDS means coding sequence (These are the only feature we're interested in)
				featQualifers = features[y].qualifiers # Each feature contains a dictionary called quailifiers which contains           
				                                       # data about the sequence feature (for example the translation)
				
				# Gets the required qualifers. Uses featQualifers.get to return the quatifer or a default value if the quatifer				# is not found. Calls strip to remove unwanted brackets and ' from quantifer before storing it as a string.
				protein_id = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')
				if protein_id == 'no_protein_id':
					continue # Skips the iteration if protien has no id.
				gene = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
				product = str(featQualifers.get('product','no_product_name')).strip('\'[]')
				translated_protein = str(featQualifers.get('translation','no_translation')).strip('\'[]')
				
				fasta.append((">" + protein_id + " " + gene + "-" + product + "\n" + translated_protein + "\n"))
	return fasta
#------------------------------------------------------------------------------------------------------------

# 4: Checks if genome is a WGSS project. 
def isSSProject(sequence):
	WGSSProjectRegex = re.compile("[a-zA-Z]{4,6}\d{8,10}")
	m = WGSSProjectRegex.match(sequence.id) 
	if m:
		matched = True
	else:
		matched = False
	return matched
#------------------------------------------------------------------------------------------------------------
# 5: When passed an array of sequence record objects returns an array of fasta strings for each annotation.
#    This implimenation is "quick and dirty" shall be replaced in later versions. 
def extractContigs(seqList):
	# Regexs for contig accession extraction (These could be pulled out of function so they are only compiled once.)
	WGSSRangeRegex = re.compile("/wgs=\['{1}[a-zA-Z0-9]{12,14}',[ ]{1}'[a-zA-Z0-9]{12,14}']{1}")
	WGSSRangeAccessionsRegex = re.compile("[a-zA-Z0-9]{12,14}")
	AccessionBaseRegex = re.compile("^[a-zA-Z]{4}\d{2}")

	SeqRecords = getSeqRecords(seqList)
	
	contigList = []

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
		
		# Creates accession list
		for x in range(WGSSRangeMin, (WGSSRangeMax + 1)):
			contigAccession = AccessionBase
			contigAccession += ("{0:0" + str(zeroOffset) + "d}").format(x) # Uses zero offset to make accessions proper length.
			contigList.append(contigAccession)
		
	return contigList