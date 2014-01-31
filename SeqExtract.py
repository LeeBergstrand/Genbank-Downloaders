#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A module that contains functions that extract and modify genbanks records. 
#
#           This module requires the Biopython module: http://biopython.org/wiki/Download
#------------------------------------------------------------------------------------------------------------
#============================================================================================================
# Imports:

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
				featQualifers = features[y].qualifiers # Each feature contains a dictionary called quailifiers which contains data about         
				                                       # the sequence feature (for example the translation)
				
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
