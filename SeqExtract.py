#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Description: A module that contains functions for extraction of genbank records.
#
# Requirements: - This module requires the Biopython module: http://biopython.org/wiki/Download
# ------------------------------------------------------------------------------------------------------------
# ============================================================================================================
# Imports and Setup:

import re
import sys

from Bio import Entrez
from Bio import SeqIO

AccessionBaseRegex = re.compile("^[a-zA-Z]{4}\d{2}")
WGSSProjectRegex = re.compile("[a-zA-Z_]{4,7}\d{8,10}")


# ============================================================================================================
# Functions:


# 1: Sets up user email.
def entrezEmail(email):
    Entrez.email = email


# ------------------------------------------------------------------------------------------------------------
# 2: When passed an array of accessions from NCBI it returns a list of sequence objects matching those accessions.
def getSeqRecords(seqList, database_type="nucleotide"):
    try:
        print("Requesting sequence data from genbank...")
        handle = Entrez.efetch(db=database_type, id=seqList, rettype="gb",
                               retmode="genbank")  # Gets records and stores them.
        print("Starting download...")
        SeqRecords = list(SeqIO.parse(handle, "genbank"))  # Creates a list of SeqRecord objects from genbank files e.
        print("Download Complete.")
        handle.close()  # Closes handle since it is no longer needed.
    except IOError:
        print("Failed to connect to NCBI server. ")
        sys.exit(1)
    return SeqRecords


# ------------------------------------------------------------------------------------------------------------
# 3: When passed a sequence record object returns an array of fasta strings for each annotation.
def getProteinAnnotationFasta(seqRecord):
    fasta = []
    features = seqRecord.features  # Each sequence has a list (called features) that stores seqFeature objects.
    for feature in features:  # For each feature on the sequence
        if feature.type == "CDS":  # CDS means coding sequence (These are the only feature we're interested in)
            featQualifers = feature.qualifiers  # Each feature contains a dictionary called qualifiers which contains
            # data about the sequence feature (for example the translation)

            # Gets the required qualifier. Uses featQualifers.get to return the qualifier or a default value if the quantifier
            # is not found. Calls strip to remove unwanted brackets and ' from qualifier before storing it as a string.
            protein_id = str(featQualifers.get('protein_id', 'no_protein_id')).strip('\'[]')
            if protein_id == 'no_protein_id':
                continue  # Skips the iteration if protein has no id.
            gene = str(featQualifers.get('gene', 'no_gene_name')).strip('\'[]')
            product = str(featQualifers.get('product', 'no_product_name')).strip('\'[]')
            translated_protein = str(featQualifers.get('translation', 'no_translation')).strip('\'[]')
            fasta.append((">" + protein_id + " " + gene + "-" + product + "\n" + translated_protein + "\n"))
    return fasta


# ------------------------------------------------------------------------------------------------------------
# 4: Passed a sequence record object returns a list of csv rows. Each row is list containing info for each annotation.
def getProteinAnnotationCSV(seqRecord):
    csvRowSet = []  # Master list of all rows.
    for feature in seqRecord.features:
        if feature.type == "CDS":
            csvRow = []  # Created new list for each row. (This is for compatibility with the csv module's row write)
            CDSLocal = feature.location  # Gets the feature location
            featQualifers = feature.qualifiers  # Gets sequence quantifiers.

            # Gets sequence quantifiers.
            gene = str(featQualifers.get('gene', 'no_gene_name')).strip('\'[]')
            product = str(featQualifers.get('product', 'no_product_name')).strip('\'[]')
            proteinID = str(featQualifers.get('protein_id', 'no_protein_id')).strip('\'[]')
            locus = str(featQualifers.get('locus_tag', 'no_locus_tag')).strip('\'[]')
            if proteinID == 'no_protein_id':
                continue  # Skips the iteration if protein has no id.
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

            csvRowSet.append(csvRow)  # Appends the csvRow list to the master list of rows.

    return csvRowSet


# ------------------------------------------------------------------------------------------------------------
# 5: Checks if genome is a WGSS project. 
def isSSProject(sequence):
    is_ssp = False

    if sequence.annotations.get('wgs'):
        is_ssp = True

    return is_ssp


# ------------------------------------------------------------------------------------------------------------
# 6: When passed a list of sequence record objects returns an list of fasta strings for each annotation.
#    This implementation is "quick and dirty" shall be replaced in later versions.
def extractContigs(seqList):
    SeqRecords = getSeqRecords(seqList)
    contigList = []

    for WGSS in SeqRecords:

        WGSSRange = WGSS.annotations["wgs"]  # Extracts WGSS contig accession range.

        if len(WGSSRange) == 2:
            # Extracts the accession base from first contig accession number.
            AccessionBase = (AccessionBaseRegex.findall(WGSSRange[0])[0])

            # Takes both the the min and max accession and slices off (using python's string slice syntax s[start:end:step])
            # the accession base code leaving the numerical difference between the contigs.
            # Converts these differences to integers.
            WGSSRangeMin = int(WGSSRange[0][6:])
            WGSSRangeMax = int(WGSSRange[1][6:])

            # WGSS accession number length actually varies. Its normally 12 characters but I have seen 13 before.
            # The code block below accounts for this.
            zeroOffset = 6
            accessionLength = len(WGSSRange[0])
            if accessionLength != 12:
                zeroOffset = accessionLength - 6  # 6 is the length of the standard accession base.

            # Creates accession list
            for x in range(WGSSRangeMin, (WGSSRangeMax + 1)):
                contigAccession = AccessionBase
                contigAccession += ("{0:0" + str(zeroOffset) + "d}").format(
                    x)  # Uses zero offset to make accessions proper length.
                contigList.append(contigAccession)
        else:
            contigList.append(WGSSRange[0])  # If one contig, simply append it to the list.

    return contigList
