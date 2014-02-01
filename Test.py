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

seqList = ["M31939"]
seqList.append("EU013931")
seqList.append("DQ345780")
#seqList.append("AVCO01000001")

entrezEmail("lee.h.bergstrand@gmail.com")
SeqRecords = getSeqRecords(seqList)

for x in SeqRecords:
	for y in x.features:
		if y.type == "CDS":
			csvRow = []
			blam = y.location
			featQualifers = y.qualifiers
			
			gene    = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
			product = str(featQualifers.get('product','no_product_name')).strip('\'[]')
			protein_id = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')

			csvRow  = [x.annotations["organism"]]
			csvRow.append(x.id)
			csvRow.append(protein_id)
			csvRow.append(gene)
			csvRow.append(str(blam.start))
			csvRow.append(str(blam.end))
			csvRow.append(str(blam.strand))
			csvRow.append(product)
			
			print csvRow