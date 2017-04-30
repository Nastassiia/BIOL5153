#!/usr/bin/env python3.6

import csv
import sys 

usage=sys.argv[0]+":genome.fasta features.gff"

if len(sys.argv)<3: 
    print(usage)
    sys.exit("\nThis script requires both a FASTA file and a GFF file\n")    

fasta_file=open(sys.argv[1], 'r')
fasta_line=fasta_file.readlines()[1]



strings_list={}
 
lines=open(sys.argv[2], 'r')
for line in lines:
    row=line.split("\t")
    key=row[2]
    start=int(row[3])-1
    end=row[4]
    if key in strings_list.keys():
        strings_list[key]=strings_list[key]+fasta_line[int(start):int(end)]
    else:
        strings_list[key]=fasta_line[int(start):int(end)]
            
       
def GC_count(string):
    count=sum(string.count(x) for x in ("G","g","C","c"))
    percent=round(count/len(string)*100,2)
    return percent
            
for key in strings_list.keys():
    print(key,'\t',len(strings_list[key]),'\t', str(round(len(strings_list[key])/len(fasta_line)*100,2))+'%'          , '\t',GC_count(strings_list[key]))
    
    
print(len(fasta_line))    


print('\n', "part of code which should extract CDS, reverse_complement '-' strands and put exons in a right order. However I did not check it thoroughly", '\n' )

import Bio
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

fasta_file=open(sys.argv[1], 'r')
fasta_line=fasta_file.readlines()[1]
cds_list={}
#here I get the first gene name with described 

lines=open(sys.argv[2], 'r')
for line in sorted(lines, key=lambda line: line.split()[10]):
    line1=line.split(" ")
    name2=line1[2] #get a gene name
    line2=line.split("\t")
    feature=line2[2]
    exons=line1[3]+line1[4]
    if (feature=='CDS' and 'exon' in exons):
        break       
        
        
exon_list={}        
lines=open("watermelon.gff", 'r')
for line in sorted(lines, key=lambda line: line.split()[10]):
    line1=line.split(" ")
    name=line1[2]
    line2=line.split("\t")
    if re.search('CDS', line):
        exons=line1[3]+line1[4] #check for multiple exons
        strand=line2[6]
        start=int(line2[3])-1 #start of strings
        end=line2[4] #end of strings
        if strand=='-':
            gene_seq=fasta_line[int(start):int(end)]
            dna=Seq(gene_seq, generic_dna)
            gene_seq=str(dna.reverse_complement())
        else:
            gene_seq=fasta_line[int(start):int(end)]
    
        if ('exon' in exons):  
            if name2==name:
                exon_list[exons]=gene_seq
                    #print(exons)
            else:
                a=sorted(exon_list.keys())
                cds_list[name2]=""
                for i in range(0,len(a)):
                    key=a[i]
                    gene_seq2=exon_list[key]
                    cds_list[name2]=(cds_list[name2]+gene_seq2)
                print(name2,'\n', cds_list[name2])        
                name2=name
                exon_list={}
                exon_list[exons]=gene_seq
                
        else:
            cds_list[name]=gene_seq
            print(name, '\n', cds_list[name])
    
