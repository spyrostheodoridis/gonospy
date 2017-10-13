import numpy as np
import os
from itertools import *

def writeLinesToFile(outfile, line_list):
    # the file should be already opened!
    for item in line_list:
        outfile.write('%s' % item)

def filterReads(indir, qThresh = 30, pQ = 0.5, trmBases = 0, repThrs = 1000):
    #lines for each read/record
    n=4
    # qThresh -> base quality threshold
    # pQ -> percentage of the read to be checked for quality score
    # trmBases -> number of bases to trim
    # repThrs -> repetitive read threshold

    #define output directory
    cleanFiles = '{}/cleanFiles'.format(indir)

    #make output directory
    os.mkdir(cleanFiles)

    #get Illumina Q-scores
    with open('{}/illumina_Q_scores_phred33.txt'.format(indir), 'r') as Qscores:
        #dictionary to store symbol and value
        illD = {}
        for line in Qscores:
            illD[line.split('\t')[0]] = line.split('\t')[1].rstrip()

    listdir = os.listdir(indir)
    nIndv = 0

    #log_file
    flog = open(os.path.join(cleanFiles,'trimmLog.txt'), 'w')
    flog.write('\t'.join(['sample', 'total_reads', 'lowQ_reads_discarded','Repetative_reads_discarded', 'Number_of_repetative_stacks', '\n']))

    for fname in listdir:
        if fname.endswith('fastq.gz'):
            os.system('gunzip -k {}/{}'.format(indir, fname))
            
            print('Processing file %s' %fname)
            fileInfo = []
            totalReads = 0
            lowQreads = 0
            repReads = 0
            repStacks = 0
            nIndv += 1
            infile = "{}/{}".format(indir, '.'.join(fname.split('.')[:-1]))
            with open(infile, 'r') as f:
                
                dic = {} #dictionary to store the id's as keys and read and quality as values
                countD = {} #dictionary to store counts (coverage) for each read
                while True:
                    nextLines = list(islice(f, n))
                    if not nextLines:
                        break

                    totalReads += 1
                    #get only the read sequence and trim it
                    rowRead = nextLines[1].rstrip()
                    read = rowRead[:-trmBases] if trmBases else rowRead
                    readLen = len(read)
                    #get the quality of the trimmed read
                    readQual = nextLines[3].rstrip()[:-trmBases] if trmBases else nextLines[3].rstrip()
                    #get read ID
                    readID = nextLines[0].rstrip()
                    
                    #mark/discard reads with at least one uncalled base. '#' is for the uncalled base
                    if '#' in readQual:
                        hiQRead = False
                    else:
                        hiQRead = True
                        
                    #use the illD dictionary to convert ascii to integer scores
                    Qint = [int(illD[i]) for i in readQual]

                    #check the quality of the first p_Q% of the reads
                    Qsub = Qint[0:(int(pQ*readLen))]
                    #check if two thirds of the bases are above the Q threshold
                    nBasesHighQ = 0
                    for score in Qsub:
                        if score >= qThresh:
                            nBasesHighQ += 1

                    #get only reads with a mean quality score above the threshold
                    if nBasesHighQ >= (2/3)*(len(Qsub)) and hiQRead == True:
                        dic[readID] = [read, readQual]
                        #count coverage
                        if read in countD:
                            #increase number of reads by 1
                            countD[read] += 1
                        else:
                            nOfReads = 1
                            countD[read] = nOfReads

                    else:
                        lowQreads += 1

                fout = open(os.path.join(cleanFiles,'{}.fastq'.format(fname.split('.')[0])), 'w')

                for key, value in dic.items(): #for each read in the dictionary
                    readStr = value[0]
                    readQuality = value[1]
                    #if read is low copy
                    if countD[readStr] <= repThrs:
                        readLines = [key, '\n', readStr, '\n', '+\n', readQual, '\n']
                        writeLinesToFile(fout, readLines)

                #get statistics for highly repetitive reads
                for read, count in countD.items():
                    if count >= repThrs:
                        repReads += count
                        repStacks += 1

                fout.close()

                logList = [str(fname.split('.')[0]),str(totalReads),str(lowQreads),str(repReads), str(repStacks), '\n']
                flog.write('\t'.join(logList))
                flog.flush()
                print('Processed Sample %s (%s) \nTotal number of reads: %s, LowQ Reads discarded: %s, Highly Repetitive Reads discarded: %s' %(fname.split('.')[0], nIndv, totalReads, lowQreads, repReads))
            os.remove(infile)
    print('Finished\n')
    flog.close()