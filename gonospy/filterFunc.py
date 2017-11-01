import numpy as np
import os
from itertools import *

def writeLinesToFile(outfile, line_list):
    # the file should be already opened!
    for item in line_list:
        outfile.write('%s' % item)

def filterReads(indir, qThresh = 30, pQ = 0.5, percQual = 2/3, trmBases = 0, repThrs = 1000):
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
                    if nBasesHighQ >= (percQual)*(len(Qsub)) and hiQRead == True:
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


def filterVCF(infile, refDepthThrs = 5, varDepthThrs = 5, maxAlleles = 1, missSampThrs = 0.2, totVarFreqThrs = 0.05, 
              totRefFreqThrs = 0.05, varCovThrs = 0.05, totCovThrs = 0.05, locus = False):
    
    # refDepthThrs -> coverage threshold for reference allele
    # varDepthThrs -> coverage threshold for the alternative allele
    # maxAlleles -> maximum number of alternative alleles
    # missSampThrs -> percentage of missing samples allowed
    # totVarFreqThrs -> minimum frequency of the variant allele across all genotypes
    # totRefFreqThrs -> minimum frequency of the reference allele across all genotypes
    # varCovThrs -> minimum frequency of the variant allele (coverage) in each genotype
    # totCovThrs -> minimum frequency of the reference allele (coverage) in each genotype
    
    FiltFile = open('filt_{}.vcf'.format(infile.split('.')[0]), 'w')

    #store the total number of missing genotypes
    totalMissData = 0
    nsn = 0
    #read file
    with open('{}'.format(infile), 'r') as VarS:
        NofSNPs = 0
        Loci = set()
        for line in VarS:
            if line.startswith('#CHROM'):
                #get samples
                header = line.strip().split('\t')
                Samples = header[9:]
            ######################################
            ########for each SNP line#############
            ######################################
            if line.startswith('un'):
                SNPLine = line.strip().split('\t')
                Locus = SNPLine[2].split('_')[0] if locus == True else SNPLine[0]
                SNP_ID = SNPLine[2]

                #Sample wide info
                NoOfVariants = len(SNPLine[4].split(','))

                if NoOfVariants > maxAlleles:
                    nsn += 1
                INFO = SNPLine[7].split(';')
                NS = [int(NS.split('=')[1]) for NS in INFO if NS.startswith('NS')][0] #Number of Samples With Data
                AF = [float(AF.split('=')[1].split(';')[0]) for AF in INFO if AF.startswith('AF')][0] #Allele Frequency of reference

                #Sample specific info            
                GenotypeList = SNPLine[9:]

                #get sample specific info
                RefSamplesNew = 0
                VarHetNew = 0
                VarHomNew = 0
                MissSampNew = 0
                GenotypeListNew = []
                #loop through all the genotypes and manipulate each genotype
                for i, GenotypeInfo in enumerate(GenotypeList):
                    #create a dictionary with the SNP position as key and the genotype as value
                    SNPdic = {}
                    #if the genotype has been called
                    if not GenotypeInfo.startswith('./.'):
                        genotype = GenotypeInfo.split(':')[0]
                        TotDepth = int(GenotypeInfo.split(':')[1]) #Read Depth
                        RefDepth = int(GenotypeInfo.split(':')[2].split(',')[0]) #Depth of reference-supporting bases
                        VarDepth = int(GenotypeInfo.split(':')[2].split(',')[1]) #Depth of variant-supporting bases
                        VarFreq = VarDepth / TotDepth #Variant allele frequency

                        #################################
                        #filter by coverage. If both states are supported by depth (heterozygous)
                        if RefDepth >= refDepthThrs and VarDepth >= varDepthThrs:

                            #################################
                            #correct for variant frequency
                            #genotype homozygous for reference
                            if VarFreq <= varCovThrs:
                                RefSamplesNew += 1
                                genotypeNew = '0/0'
                                GenotypeInfoNew = GenotypeInfo.replace(genotype, genotypeNew)
                                #replace genotype to original
                                GenotypeListNew.append(GenotypeList[i].replace(GenotypeInfo,GenotypeInfoNew))
                                SNPdic[SNP_ID] = [genotypeNew, TotDepth, RefDepth, VarDepth, VarFreq]
                            #genotype homozygous for variant
                            elif VarFreq >= 1 - varCovThrs:
                                VarHomNew += 1
                                hetState = genotype.split('/')[1]
                                genotypeNew = '%s/%s' % (hetState, hetState)
                                GenotypeInfoNew = GenotypeInfo.replace(genotype, genotypeNew)
                                #replace genotype to original
                                GenotypeListNew.append(GenotypeList[i].replace(GenotypeInfo,GenotypeInfoNew))
                                SNPdic[SNP_ID] = [genotypeNew, TotDepth, RefDepth, VarDepth, VarFreq]
                            #genotype heterozygous
                            elif varCovThrs < VarFreq < 1 - varCovThrs:
                                VarHetNew += 1
                                RefState = 0
                                hetState = genotype.split('/')[1]
                                genotypeNew = '%s/%s' % (RefState, hetState)
                                GenotypeInfoNew = GenotypeInfo.replace(genotype, genotypeNew)
                                #replace genotype to original
                                GenotypeListNew.append(GenotypeList[i].replace(GenotypeInfo,GenotypeInfoNew))
                                SNPdic[SNP_ID] = [genotype, TotDepth, RefDepth, VarDepth, VarFreq]
                            else:
                                print('opps, check {} SNP for errors'.format(SNP_ID))

                        #variant coverage lower than threshold (homozygous for reference)
                        elif RefDepth >= refDepthThrs and VarDepth < varDepthThrs:
                            RefSamplesNew += 1
                            genotypeNew = '0/0'
                            GenotypeInfoNew = GenotypeInfo.replace(genotype, genotypeNew)
                            #replace genotype to original
                            GenotypeListNew.append(GenotypeList[i].replace(GenotypeInfo,GenotypeInfoNew))
                            SNPdic[SNP_ID] = [genotypeNew, TotDepth, RefDepth, VarDepth, VarFreq]

                        #reference coverage lower than threshold
                        elif RefDepth < refDepthThrs and VarDepth >= varDepthThrs:
                            VarHomNew += 1
                            hetState = genotype.split('/')[1]
                            genotypeNew = '%s/%s' % (hetState, hetState)
                            GenotypeInfoNew = GenotypeInfo.replace(genotype, genotypeNew)
                            #replace genotype to original
                            GenotypeListNew.append(GenotypeList[i].replace(GenotypeInfo,GenotypeInfoNew))
                            SNPdic[SNP_ID] = [genotypeNew, TotDepth, RefDepth, VarDepth, VarFreq]

                        #both states have lower coverage that the theshold    
                        else:
                            MissSampNew += 1
                            genotypeNew = './.'
                            GenotypeInfoNew = GenotypeInfo.replace(genotype, genotypeNew)
                            #replace genotype to original
                            GenotypeListNew.append(GenotypeList[i].replace(GenotypeInfo,GenotypeInfoNew))
                            SNPdic[SNP_ID] = [genotypeNew, TotDepth, RefDepth, VarDepth, VarFreq]

                    else: #missing genotype
                        GenotypeListNew.append(GenotypeList[i])
                        MissSampNew += 1
                        SNPdic[SNP_ID] = 'missing'

                #the new frequency of the variant. If none of the genotypes passed the coverage thresholds then  TotVarFreq = 0
                try:
                    TotVarFreq = (VarHetNew + (VarHomNew*2)) / ((RefSamplesNew + VarHetNew + VarHomNew)*2)
                except:
                    TotVarFreq = 0

                #the new frequency of the reference. If none of the genotypes passed the coverage thresholds then TotRefFreq = 0
                try:
                    TotRefFreq = (VarHetNew + (RefSamplesNew*2)) / ((RefSamplesNew + VarHetNew + VarHomNew)*2)
                except:
                    TotRefFreq = 0
                ##############################################
                if MissSampNew/len(Samples) <= missSampThrs and TotVarFreq >= totVarFreqThrs and TotRefFreq >= totRefFreqThrs and NoOfVariants <= maxAlleles:
                    #new info based on the replacement of the genotypes
                    infoNew = ['NS={:d}'.format(len(Samples) - MissSampNew), 'AF={:.3f}'.format(TotVarFreq)]
                    #create new line
                    lineNew = line.replace(SNPLine[7], ';'.join(infoNew)).replace('\t'.join(GenotypeList), '\t'.join(GenotypeListNew))
                    #write line that has passed all above filters (TotVarFreq, Miss) to the new file
                    FiltFile.write(lineNew)
                    NofSNPs += 1
                    Loci.add(Locus)

                    totalMissData += MissSampNew
            ######################################
            ######################################
            ######################################
            #write the first lines of the vcf to the new file
            else:
                FiltFile.write(line)

            #clean memory
            try:
                del(line, GenotypeList, GenotypeListNew)
                gc.collect()
            except:
                continue

    FiltFile.close()
    print('Number of SNPs: %s' %NofSNPs)
    print('Number of Loci / Chromosomes: %s' %len(Loci))
    print('Total Percentage of Missing Data: %.2f' %(totalMissData/(len(Samples)*NofSNPs)))

