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

seqList = ["CBMO000000000"]
seqList.append("BAEC00000000")

entrezEmail("leemacboy@gmail.com")

print extractContigs(seqList)




	

