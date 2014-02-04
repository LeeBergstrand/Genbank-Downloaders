#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: 
#
#           This script requires the Biopython module: http://biopython.org/wiki/Download
#  
# Usage: 
# Example:
#----------------------------------------------------------------------------------------

from SeqExtract import entrezEmail	
from SeqExtract import extractContigs
from SeqExtract import getSeqRecords
from SeqExtract import isSSProject
from SeqExtract import getProtienAnnotationFasta

def getProtienAnnotationCSV(seqRecord):	
	csvRowSet = [] # Master list of all rows.
	for feature in seqRecord.features:
		if feature.type == "CDS":
			csvRow = [] # Created new list for each row. (This is for compatibility with the csv module's row write)
			CDSLocal = feature.location # Gets the feature location
			featQualifers = feature.qualifiers # Gets sequence quantifiers.
			
			# Gets sequence quantifiers.
			gene       = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
			product    = str(featQualifers.get('product','no_product_name')).strip('\'[]')
			proteinID  = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')
			locus      = str(featQualifers.get('locus_tag','no_locus_tag')).strip('\'[]')
			if proteinID == 'no_protein_id':
				continue # Skips the iteration if protien has no id.
			# Append quantifers and other information to the csv row list.
			csvRow.append(seqRecord.annotations["organism"])
			csvRow.append(seqRecord.id)
			csvRow.append(proteinID)
			csvRow.append(gene)
			csvRow.append(str(CDSLocal.start))
			csvRow.append(str(CDSLocal.end))
			csvRow.append(str(CDSLocal.strand))
			csvRow.append(locus)
			csvRow.append(product)
			
			csvRowSet.append(csvRow) # Appends the csvRow list to the master list of rows.
	
	return csvRowSet

seqList = ["CP000431"]

entrezEmail("")
SeqRecords = getSeqRecords(seqList)

for sequence in SeqRecords:
	csvRows = getProtienAnnotationCSV(sequence)
	for x in csvRows:
		print x

# def getProtienAnnotationFasta(seqRecord):
# fasta = []
# features = seqRecord.features # Each sequence has a list (called features) that stores seqFeature objects.
# for feature in features: # For each feature on the sequence
#	if feature.type == "CDS": # CDS means coding sequence (These are the only feature we're interested in)
#		featQualifers = feature.qualifiers # Each feature contains a dictionary called quailifiers which contains           
#		                                   # data about the sequence feature (for example the translation)
#		
#		# Gets the required qualifers. Uses featQualifers.get to return the quatifer or a default value if the quatifer			# is not found. Calls strip to remove unwanted brackets and ' from quantifer before storing it as a string.
#		protein_id = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')
#		if protein_id == 'no_protein_id':
#			continue # Skips the iteration if protien has no id.
#		gene    = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
#		product = str(featQualifers.get('product','no_product_name')).strip('\'[]')
#		translated_protein = str(featQualifers.get('translation','no_translation')).strip('\'[]')
#		fasta.append((">" + protein_id + " " + gene + "-" + product + "\n" + translated_protein + "\n"))
# return fasta