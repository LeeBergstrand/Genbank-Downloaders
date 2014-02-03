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
from SeqExtract import getProtienAnnotationFasta2

seqList = ["M31939"]
seqList.append("EU013931")
seqList.append("DQ345780")
seqList.append("AVCO01000001")

entrezEmail("")
SeqRecords = getSeqRecords(seqList)

print getProtienAnnotationFasta()