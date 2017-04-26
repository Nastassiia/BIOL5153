#!/usr/bin/env python3.6

import csv
import sys 

usage=sys.argv[0]+":genome.fasta features.gff"

if len(sys.argv)<3: 
    print(usage)
    sys.exit("\nThis script requires both a FASTA file and a GFF file\n")
    

fasta_file=open(sys.argv[1], 'r')
fasta_line=fasta_file.readlines()[1]

names_list=['CDS', 'intron', 'misc_feature', 'rRNA', 'repeat_region', 'tRNA']
strings_list=['exon_string', 'intron_string', 'misc_string','rRNA_string', 'repeats_string', 'tRNA_string']
for i in range(0,len(strings_list)):
    strings_list[i]=""

with open(sys.argv[2], 'r') as csvfile:
    line=csv.reader(csvfile, delimiter="\t")
    for row in line:
        if row[2]==names_list[0]:
            #print(row)
            start=int(row[3])-1
            end=row[4]
            strings_list[0]=strings_list[0]+fasta_line[int(start):int(end)]
            
        if row[2]==names_list[1]:
            start=int(row[3])-1
            end=row[4]
            strings_list[1]=strings_list[1]+fasta_line[int(start):int(end)]
            
        if row[2]==names_list[2]:
            start=int(row[3])-1
            end=row[4]
            strings_list[2]=strings_list[2]+fasta_line[int(start):int(end)]
        if row[2]==names_list[3]:
            start=int(row[3])-1
            end=row[4]
            strings_list[3]=strings_list[3]+fasta_line[int(start):int(end)] 
        if row[2]==names_list[4]:
            start=int(row[3])-1
            end=row[4]
            strings_list[4]=strings_list[4]+fasta_line[int(start):int(end)]
        if row[2]==names_list[5]:
            start=int(row[3])-1
            end=row[4]
            strings_list[5]=strings_list[5]+fasta_line[int(start):int(end)]   
def GC_count(string):
    count=sum(string.count(x) for x in ("G","g","C","c"))
    percent=round(count/len(string)*100,2)
    return percent


def N_count(string, N):
    l=len(string)
    string=string.upper()
    numb=round(string.count(N)/l,2)
    return numb
    
            
for i in range(0,len(strings_list)):
    print(names_list[i],'\t',len(strings_list[i]),'\t', str(round(len(strings_list[i])/len(fasta_line)*100,2))+'%' , '\t',GC_count(strings_list[i]))
    print(names_list[i], '\t', "A content:", N_count(strings_list[i], "A"), '\t',  "T content:", N_count(strings_list[i], "T"), '\t',  "G content:", \
N_count(strings_list[i], "G"), '\t',  "C content:", N_count(strings_list[i], "C") )

    
print(len(fasta_line))    

            
        
        




