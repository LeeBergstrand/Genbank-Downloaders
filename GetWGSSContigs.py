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

seqList = ["CBMO000000000"]
seqList.append("BAEC00000000")
seqList.append("CP003080")

entrezEmail("")
SeqRecords = getSeqRecords(seqList)

for x in SeqRecords:
	print x.id
	print isSSProject(x)
#print extractContigs(seqList)






	

