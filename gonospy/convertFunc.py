import random
import sys
import re

# two functions for natural sorting
def strToInt(text):
    return int(text) if text.isdigit() else text

def naturalOrder(text):
    return [ strToInt(c) for c in re.split('(\d+)', text) ]

def converVCF(infile, outFormat, missChr = '?', oneSNP = False, nOfRep = 1, excludeInv = False, randomAllele = False, locusID = True, populationFile = None):

    #dictionary to store SNPs per individual
    IndvSNPs = {}

    with open(infile, 'r') as VarS:

        #number of sites
        nSites = 0
        #number of loci
        Loci = {}
        #SNP names
        allSNPs = []

        for line in VarS:
            if line.startswith('#CHROM'):
                #get samples
                header = line.strip().split('\t')
                Samples = header[9:]
                #create sample:population dictionary
                if populationFile != None:
                    with open(populationFile, 'r') as pF:
                        samPopD = {}
                        for line in pF:
                            samPopD[line.strip().split('\t')[0]] = line.strip().split('\t')[1]
                        Pop_Index = { pop:i for i, pop in enumerate(sorted(set(samPopD.values()), key = naturalOrder))}
                        PopIndex = { samp:Pop_Index[pop] for samp, pop in samPopD.items() }
                        pF.close()
                else: #get population names from the sample names
                    Populations = set()
                    for sample in Samples:
                        popName = sample.split('_')[-1]
                        Populations.add(popName)
                    Populations = sorted(Populations, key = naturalOrder)
                    Pop_Index = { pop:i for i, pop in enumerate(Populations) }
                    PopIndex = { samp:Pop_Index[samp.split('_')[-1]] for samp in Samples }
            #get genotypes
            if line.startswith('un'):
                nSites += 1

                SNPLine = line.strip().split('\t')
                ref = SNPLine[3]
                alt = SNPLine[4]
                #Variant-wide genotype info
                GenotypeL = [sample.split(':')[0] for sample in SNPLine if sample.startswith(('./', '0/', '1/'))] # or the max alleles per SNP
                GenotypeList = [w.replace('.', missChr) for w in GenotypeL]
                GenotypeList = [w.replace('0', ref) for w in GenotypeList]
                GenotypeList = [w.replace('1', alt) for w in GenotypeList]

                #check if snp is invariant among genotypes
                genNoMiss = [i for i in GenotypeList if missChr not in i]
                distinctGenotypes = list(set(genNoMiss))
                if len(distinctGenotypes) < 3:
                    if distinctGenotypes[0].split('/')[0] in distinctGenotypes[1] or distinctGenotypes[0].split('/')[1] in distinctGenotypes[1]:
                        siteGen = 'invariant'
                    else:
                        siteGen = 1
                else:
                    siteGen = 1

                if excludeInv == True and siteGen != 'invariant' or excludeInv == False:
                    if locusID == True:
                        Locus = SNPLine[2].split('_')[0]
                    else:
                        Locus = SNPLine[0]
                    SNP_ID = SNPLine[2]
                    allSNPs.append(SNP_ID)
                    Loci.setdefault(Locus, []).append(SNP_ID)

                    for i, genot in enumerate(GenotypeList):
                        AllelesPerIndv = genot.split('/')
                        sampl = Samples[i]
                        IndvSNPs.setdefault(sampl, {}).setdefault(Locus, {})[SNP_ID] = AllelesPerIndv

    VarS.close()                

    if oneSNP == True: #randomly choose one SNP per locus

        for rep in range(nOfRep):
            if outFormat == 'fastStructure':
                outfile = open('{}_{}.fastStructure.str'.format(infile.split('.')[:-1][0],rep+1), 'w')
            elif outFormat == 'structure':
                outfile = open('{}_{}.structure'.format(infile.split('.')[:-1][0], rep+1), 'w')
            elif outFormat == 'phylip':
                outfile = open('{}_{}.phy'.format(infile.split('.')[:-1][0], rep+1), 'w')

            #pick one random SNP per locus
            RandSNP = []
            for locus in sorted(Loci, key = naturalOrder):
                random_SNP = random.choice(Loci[locus])
                RandSNP.append(random_SNP)

            if outFormat == 'structure':
                RandSNP.insert(0, '\t')
                outfile.write('\t'.join(RandSNP))
                outfile.write('\n')
            if outFormat == 'phylip':
                outfile.write('\t{}\t{}\n'.format(len(IndvSNPs), len(RandSNP)))

            samp = 0
            for Indv in sorted(IndvSNPs, key = naturalOrder):
                samp += 1            
                sys.stdout.write('\rWriting SNPs for individual {} {}/{} rep {}'.format(Indv, samp, len(IndvSNPs), rep))
                pop = str(PopIndex[Indv])

                if outFormat == 'fastStructure':
                    outfile.write('\t'.join([Indv,pop,'NA\t'*4,]))
                elif outFormat == 'structure':
                    outfile.write('\t'.join([Indv,pop,'']))
                elif outFormat == 'phylip':
                    outfile.write('{}\t'.format(Indv))

                if outFormat == 'structure' or outFormat == 'fastStructure':

                    for c, CLocus in enumerate(sorted(IndvSNPs[Indv], key = naturalOrder)):
                        for snp in sorted(IndvSNPs[Indv][CLocus], key = naturalOrder):
                            if snp in RandSNP:
                                outfile.write(IndvSNPs[Indv][CLocus][snp][0])
                                if c < len(IndvSNPs[Indv])-1:
                                    outfile.write('\t')
                    outfile.write('\n')

                    if outFormat == 'fastStructure':
                        outfile.write('\t'.join([Indv,pop,'NA\t'*4,]))
                    else:
                        outfile.write('\t'.join([Indv,pop,'']))
                    #second line for each individual   
                    for c, CLocus in enumerate(sorted(IndvSNPs[Indv], key = naturalOrder)):
                        for snp in sorted(IndvSNPs[Indv][CLocus], key = naturalOrder):
                            if snp in RandSNP:
                                outfile.write(IndvSNPs[Indv][CLocus][snp][1])
                                if c < len(IndvSNPs[Indv])-1:
                                    outfile.write('\t')
                    outfile.write('\n')

                elif outFormat == 'phylip':
                    for c, CLocus in enumerate(sorted(IndvSNPs[Indv], key = naturalOrder)):
                        for snp in sorted(IndvSNPs[Indv][CLocus], key = naturalOrder):
                            if snp in RandSNP:
                                if randomAllele == True:
                                    outfile.write(IndvSNPs[Indv][CLocus][snp][random.randint(0,1)])
                                elif randomAllele == False:
                                    site = IndvSNPs[Indv][CLocus][snp]
                                    if len(set(site)) > 1:
                                        if 'A' in site and 'G' in site:
                                            outfile.write('R')
                                        if 'C' in site and 'T' in site:
                                            outfile.write('Y')
                                        if 'G' in site and 'C' in site:
                                            outfile.write('S')
                                        if 'A' in site and 'T' in site:
                                            outfile.write('W')
                                        if 'G' in site and 'T' in site:
                                            outfile.write('K')
                                        if 'A' in site and 'C' in site:
                                            outfile.write('M')
                                    else:
                                        outfile.write(IndvSNPs[Indv][CLocus][snp][0])                          
                    outfile.write('\n')

        outfile.close()


    if oneSNP == False:

        if outFormat == 'fastStructure':
                outfile = open('{}_allSNPs.fastStructure.str'.format(infile.split('.')[:-1][0]), 'w')
        elif outFormat == 'structure':
            outfile = open('{}_allSNPs.structure'.format(infile.split('.')[:-1][0]), 'w')
        elif outFormat == 'phylip':
            outfile = open('{}_allSNPs.phy'.format(infile.split('.')[:-1][0]), 'w')

        allSNPs = sorted(allSNPs, key = naturalOrder)

        if outFormat == 'structure':
            allSNPs.insert(0, '\t')
            allSNPs.append('\n')
            outfile.write('\t'.join(allSNPs))
        if outFormat == 'phylip':
            outfile.write('\t{}\t{}\n'.format(len(IndvSNPs), len(allSNPs)))

        samp = 0
        for Indv in sorted(IndvSNPs, key = naturalOrder):
            samp += 1

            sys.stdout.write('\rWriting SNPs for individual {} {}/{}'.format(Indv, samp, len(IndvSNPs)))

            pop = str(PopIndex[Indv])

            if outFormat == 'fastStructure':
                    outfile.write('\t'.join([Indv,pop,'NA\t'*4,]))
            elif outFormat == 'structure':
                outfile.write('\t'.join([Indv,pop,'']))
            elif outFormat == 'phylip':
                outfile.write('{}\t'.format(Indv))

            if outFormat == 'structure' or outFormat == 'fastStructure':

                for CLocus in sorted(IndvSNPs[Indv], key = naturalOrder):
                    for snp in sorted(IndvSNPs[Indv][CLocus], key = naturalOrder):
                        outfile.write(IndvSNPs[Indv][CLocus][snp][0])
                        outfile.write('\t')
                outfile.write('\n')
                if outFormat == 'fastStructure':
                    outfile.write('\t'.join([Indv,pop,'NA\t'*4,]))
                else:
                    outfile.write('\t'.join([Indv,pop,'']))

                #second line for each individual
                for CLocus in sorted(IndvSNPs[Indv], key = naturalOrder):
                    for snp in sorted(IndvSNPs[Indv][CLocus], key = naturalOrder):
                        outfile.write(IndvSNPs[Indv][CLocus][snp][1])
                        outfile.write('\t')
                outfile.write('\n')

            elif outFormat == 'phylip':
                for CLocus in sorted(IndvSNPs[Indv], key = naturalOrder):
                    for snp in sorted(IndvSNPs[Indv][CLocus], key = naturalOrder):
                        if randomAllele == True:
                            outfile.write(IndvSNPs[Indv][CLocus][snp][random.randint(0,1)])
                        elif randomAllele == False:
                            site = IndvSNPs[Indv][CLocus][snp]
                            if len(set(site)) > 1:
                                if 'A' in site and 'G' in site:
                                    outfile.write('R')
                                if 'C' in site and 'T' in site:
                                    outfile.write('Y')
                                if 'G' in site and 'C' in site:
                                    outfile.write('S')
                                if 'A' in site and 'T' in site:
                                    outfile.write('W')
                                if 'G' in site and 'T' in site:
                                    outfile.write('K')
                                if 'A' in site and 'C' in site:
                                    outfile.write('M')
                            else:
                                outfile.write(IndvSNPs[Indv][CLocus][snp][0])     
                outfile.write('\n')

    outfile.close()

    del GenotypeList, allSNPs, IndvSNPs
    print('\n{} file created! Random SNP: {}; Exclude Invariant Sites: {}'.format(outFormat, oneSNP, excludeInv))