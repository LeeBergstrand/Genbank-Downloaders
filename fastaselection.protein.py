#usage
#python fastaselection2.py <originalfile.fasta> <namelist>
# 	0			1		2		

import sys, screed

#inputs
filein=open(sys.argv[1],'r')
filelist=open(sys.argv[2],'r')

#outputs
outy=sys.argv[1]
out1=outy+'.subselection'
fileout=open(out1,'w')


#create a list with the names of the sequences requested
requestedsequences=[]
for line in filelist:
   line=line.strip('\n').strip('\r')
   requestedsequences.append(line)

number_records=len(requestedsequences)
print "%s records requested" % number_records


#read file, read each record, if name is in list write it, otherwise continue
counter=1
found_counter=0
found_list=[]
for record in screed.open(sys.argv[1]):
   sequence_name=record.name			#get sequence name
   if sequence_name in requestedsequences:
      found_list.append(sequence_name)
      found_counter=found_counter+1
      print "%s of %s records found" %(counter, number_records)
      sequence=record.sequence
      sequence=sequence.strip('*') #get rid of stop codon marked as *
      description=record.description
      fileout.write(">%s %s\n%s\n" %(sequence_name, description, sequence))
      counter=counter+1
   else:
      continue

print "\n"
print "%s records requested" % number_records
print "%s records found" % found_counter

if number_records == found_counter:
   print 'All fine'
else:
   print "Not found:"
   for member in requestedsequences:
      if member in found_list:
         continue
      else:
         print member
fileout.close()
filein.close()
filelist.close()
