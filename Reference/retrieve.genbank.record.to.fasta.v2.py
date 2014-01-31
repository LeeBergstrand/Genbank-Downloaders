#usage 
#python retrieve.genbank.record.py <list.of.sequences.requested> <output.fa>
# 	0					1			2

import sys
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "carden24@mail.ubc.ca"

#inputs
filelist=open(sys.argv[1],'r')

#output
fileout=open(sys.argv[2],'w')


#create a list with the names of the sequences requested
requestedsequences=[]
for line in filelist:
   line=line.strip('\n')
   requestedsequences.append(line)
print "%d Sequences requested" % len(requestedsequences)
print requestedsequences


handle = Entrez.efetch(db="protein", id=requestedsequences, rettype="gb", retmode="genbank")
records=SeqIO.parse(handle,"genbank")

for record in records:
   feat=record.features
   for f in feat:
      if f.type=="CDS":
         quali=f.qualifiers
         gene=str(quali.get('gene','no_gene_name'))
         gene=gene.strip('\'[]')
         print gene
         product=str(quali.get('product','no_product_name'))
         product=product.strip('\'[]')
#         print product
#         description=gene+'-'+product
         protein_id=str(quali.get('protein_id','no_protein_id'))
#         protein_id=str(f.qualifiers['protein_id'])
         protein_id=protein_id.strip('\'[]')
         translated_protein=str(quali.get('translation','no_translation'))
#         translated_protein=str(f.qualifiers['translation'])
         translated_protein=translated_protein.strip('\'[]')
         if protein_id=='no_protein_id':
            continue
         else:
            fileout.write(">%s %s\n%s\n" %(protein_id,gene,translated_protein))
      else:
         continue
fileout.close()

