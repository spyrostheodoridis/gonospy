import numpy as np
import os
from itertools import *
import time

start = time.time()

def write_lines_to_file(outfile, line_list):
    # the file should be already opened!
    for item in line_list:
        outfile.write("%s" % item)

#lines for each read/record
n=4

#Base quality threshold
Q_thresh = 30

#percentage of the read to be checked for quality score
p_Q = 0.5

#number of bases to trim
trm_bases = 11

#repetitive read threshold
rep_Thrs = 1000

#define input directory
indir = '/Volumes/3.6_TB_LaCie/study_2/seq/test'

#define output directory
clean_files = '/Volumes/3.6_TB_LaCie/study_2/seq/test/clean_files'

#make output directory
os.mkdir(clean_files)

#get Illumina Q-scores
with open("./illumina_Q_scores_phred33.txt", "r") as Qscores:
    #dictionary to store symbol and value
    ill_d = {}
    for line in Qscores:
        ill_d[line.split("\t")[0]] = line.split("\t")[1].rstrip()
        
#get population codes for individuals
with open("./popmap_codes.txt", "r") as PopCodes:
    Pop_D = {}
    for line in PopCodes:
        Pop_D[line.split("\t")[0]] = line.split("\t")[1]

listdir = os.listdir(indir)
n_indv = 0

#log_file
flog = open(os.path.join(clean_files,'trimm_log.txt'), 'w')
flog.write('\t'.join(['sample', 'total_reads', 'lowQ_reads_discarded','Repetative_reads_discarded', 'Number_of_repetative_stacks', '\n']))

for fname in listdir:
    if fname.startswith('pf') and fname.endswith("fastq"):
        #os.system('sudo gunzip -k %s' %fname)
        print('Processing file %s' %fname)
        file_info = []
        total_reads = 0
        low_Q_reads = 0
        repReads = 0
        repStacks = 0
        n_indv += 1
        
        extrFile = fname#[:-3]
        with open(os.path.join(indir, extrFile), 'r') as f:
            dic = {} #dictionary to store the id's as keys and read and quality as values
            countD = {} #dictionary to store counts (coverage) for each read
            while True:
                next_n_lines = list(islice(f, n))
                if not next_n_lines:
                    break
                    
                total_reads += 1
                #get only the read sequence and trim it
                row_read = next_n_lines[1].rstrip()
                read = row_read[:-trm_bases]
                readLen = len(read)
                #get the quality of the trimmed read
                read_Qual = next_n_lines[3].rstrip()[:-trm_bases]
                #get read ID
                read_ID = next_n_lines[0].rstrip()
                
                #mark/discard reads with at least one uncalled base. '#' is for the uncalled base
                if '#' in read_Qual:
                    hiQRead = False
                else:
                    hiQRead = True
                
                #use the ill_d dictionary to convert ascii to integer scores
                Q_int = [int(ill_d[i]) for i in read_Qual]
                
                #check the quality of the first p_Q% of the reads
                Q_sub = Q_int[0:(int(p_Q*readLen))]
                #check if two thirds of the bases are above the Q threshold
                n_bases_highQ = 0
                for score in Q_sub:
                    if score >= Q_thresh:
                        n_bases_highQ += 1
                
                #get only reads with a mean quality score above the threshold
                if n_bases_highQ >= (2/3)*(len(Q_sub)) and hiQRead == True:
                    dic[read_ID] = [read, read_Qual]
                    #count coverage
                    if read in countD:
                        #increase number of reads by 1
                        countD[read] += 1
                    else:
                        nOfReads = 1
                        countD[read] = nOfReads
                        
                else:
                    low_Q_reads += 1
                    
            #open file to write filtered data
            popCode = Pop_D[extrFile[0:-6]]
            fullName = '_'.join([extrFile[0:-6], popCode])
            fout = open(os.path.join(clean_files,'%s.fastq' %fullName), 'w')
            
            for key, value in dic.items(): #for each read in the dictionary
                readStr = value[0]
                read_Quality = value[1]
                #if read is low copy
                if countD[readStr] <= rep_Thrs:
                    read_lines = [key, '\n', readStr, '\n', '+\n', read_Quality, '\n']
                    write_lines_to_file(fout, read_lines)
                
            #get statistics for highly repetitive reads
            for read, count in countD.items():
                if count >= rep_Thrs:
                    repReads += count
                    repStacks += 1
            
            fout.close()
            
            log_list = [str(fullName),str(total_reads),str(low_Q_reads),str(repReads), str(repStacks), '\n']
            flog.write('\t'.join(log_list))
            flog.flush()
            print('Processed Sample %s (%s) \nTotal number of reads: %s, LowQ Reads discarded: %s, Highly Repetitive Reads discarded: %s' %(fullName, n_indv, total_reads, low_Q_reads, repReads))
            os.system('sudo rm %s' %extrFile)
print('Finished\n')
flog.close()