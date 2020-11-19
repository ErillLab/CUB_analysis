# -*- coding: utf-8 -*-
# Copyright Ivan Erill (erill@umbc.edu)

"""
CUB analyzer

Analyzes a set of files using the nRCA and CAI CUB indices
rRCA and CAI functions from [CodonUsageIndiceslib]

Generates, for each file in input folder, a CSV file in output folder
containing the extracted genome features and the nRCA/CAI values.

Files in input folder should be GenBank files GBK

Input: CUB_config.json JSON file containing
    - Reference set file name [a FASTA file with multiple CDS records]
    - Input folder path [a folder containing GenBank formatted genome files]
    - Output folder path [a folder were CSV files will be placed]
"""

import CodonUsagelib
import os
import sys
import json

from Bio import SeqIO

""" 
    **************************************************************************
    Auxiliary functions
    Defines auxiliary functions for processing of genome records
    **************************************************************************
"""
def readGBK(filename):
    """Reads a GBK file and returns the object handle"""
    gnome_record = SeqIO.read(filename, "genbank")
    return gnome_record

def GBK_CUB_analysis(genomeobject, CAIind, nRCAind, outfilename):
    """Takes in a GBK record object, a codon usage index and an output file.
       It goes through every CDS in the genome file, and scores its codon usage
       according to the provided index.
       
       Inputs:
       - genomeobject
           a parsed GBK genome object, as done by SeqIO.read (e.g. readGBK)
       - CAIind
           instantiated CAI index
       - nRCAind
           instantiated nRCA index
       - outfilename
           the output file name
    """

    with open(outfilename,"w") as out_handle:
        out_handle.write('GenomeID,LocusTag,ProtID,Location,CAI,nRCA,Product')
        out_handle.write('\n')
        ft_cnt=1
        #iterate features
        for feat in genomeobject.features:
            #get type features
            if feat.type == 'CDS':
                cai_val = CAIind.cai_for_gene(str(feat.location.extract(genomeobject).seq.upper()))
                nrca_val = nrca_index.nrca_for_gene(str(feat.location.extract(genomeobject).seq.upper()))

                out_handle.write(genomeobject.id + ',')
                if 'locus_tag' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['locus_tag'][0] + ',')
                elif 'note' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['note'][0] + ',')
                else:
                    out_handle.write('NO_LOCUS_TAG_' + str(ft_cnt) + ',')
                if 'protein_id' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['protein_id'][0] + ',')
                else:
                    out_handle.write(',')
                out_handle.write(str(feat.location.nofuzzy_start) + '-')
                out_handle.write(str(feat.location.nofuzzy_end) + ' : ')
                out_handle.write(str(feat.location.strand))
                out_handle.write(',')
                out_handle.write(str(cai_val))
                out_handle.write(',')
                out_handle.write(str(nrca_val))
                out_handle.write(',')
                if 'product' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['product'][0])
                out_handle.write('\n')
                ft_cnt=ft_cnt+1
            
    return(0)


#read configuration from CUB_config.json
with open('CUB_config.json') as json_file:
    settings = json.load(json_file)
    
"""Instantiate CUB objects with provided reference files
"""
#make nRCA and CAI index objects
nrca_index = CodonUsagelib.normRelativeCodonAdaptationIndex()
cai_index = CodonUsagelib.CodonAdaptationIndex()

#generate an index from reference set file specified in CUB_config.json
if os.path.exists(settings["refset_file"]):
    nrca_index.generate_index(settings["refset_file"])
    cai_index.generate_index(settings["refset_file"])
else:
    print("Cannot find the reference set file \n \
          Make sure it is in the appropriate folder")
    sys.exit()

#for each genbank file in input folder, call GCK_CUB to get and save
#CAI and nRCA index values for all CDS features
if os.path.exists(settings["input_path"]):
    dir_contents=os.listdir(settings["input_path"])
    for genomefile in dir_contents:
        print "Processing: ", genomefile
        try:
            genome = readGBK(settings["input_path"]+'/'+genomefile)
        except:
            print "Error processing: ", genomefile
            print "Revise file format. Exiting..."
            sys.exit()
        GBK_CUB_analysis(genome, cai_index, nrca_index, \
                         settings["output_path"]+'/'+ \
                         genomefile.split('.')[0]+'.csv')
else:
    print("Please make sure that the specified input folder exists")
    sys.exit()

