#! /usr/bin/env python3
# written by Noelle Schlenk (n752s925@ku.edu)
#   Filter.MQADPAR written by C. Wessinger and J. Kelly, modified by N. Schlenk
#   Filter.BestSNP written by J. Kelly, modified by N. Schlenk

import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np

class FilterVCF: # no dependencies 
    def __init__(self, out_dir = None, nF2s = 238, scaffold_index = 12):        
        # out_dir = directory for outfiles (if None, outfiles are written to working directory)
        # nF2s = number of F2 individuals in the VCF
        # scaffold_index = index (base 0) within the title of each scaffold which can be used to define that scaffold
        #   this is used to name files with a single digit rather than the entire scaffold name
        #   for example, my scaffold 0 title starts with 'PGA_scaffold0...', so I only use the character in position 12 (0) to define it             
        if out_dir != None:
            try:                                                                                    # generate directory (does not overwrite existing directory)
                self.out_dir = str(out_dir)
                os.makedirs(self.out_dir)       
            except FileExistsError:
                pass
        else:
            self.out_dir = '.'
        try:                                                                                        # verify that nF2s is an integer
            self.nF2s = int(nF2s)        
        except TypeError:
            print('NF2s must be an integer.')
            quit()            
        try:                                                                                        # verify that scaffold_index is an integer
            self.scaffold_index = int(scaffold_index)
        except TypeError:
            print('Scaffold index must be an integer.')
            quit()
    def MQADPAR(self, infile_path, outfile_path, minIndiv = 50, minMQ = 30):  # filter for MQ score, minIndiv, and informative parents      
        # written by CAW and JKK, modified by NCS
        # minIndiv = minimum number of individuals carrying SNP required to pass filter
        # minMQ = minimum MQ score required for a SNP to pass filter
        if os.path.isfile(self.out_dir + '/' + outfile_path) == True:                               # check if the outfile exists, give warning
            a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
            if a == 'n':
                quit()
        outfile = open(self.out_dir + '/' + outfile_path, 'w')
        infile = open(infile_path, 'r')
        for line in infile:                                                                         # extract MQ score and alt_allele from each line, keep line if it passes minimum                
            try:
                mq_score = [float(val.split('=')[1]) for val in line.split('\t')[7].split(';') if val.split('=')[0] == 'MQ'][0]
                alt_allele = line.split('\t')[4].split(',')
                if mq_score >= minMQ and len(alt_allele) == 1:
                    calls = 0
                    for i in range(9, 11 + self.nF2s):
                        if len(line.split('\t')[i].split(':')) < 3:
                            pass
                        else:
                            calls += 1                                                              # count number of calls
                    p1= line.replace('\n','').split('\t')[self.nF2s + 9].split(':')
                    p2 = line.replace('\n','').split('\t')[self.nF2s + 10].split(':')               # p1 must be 1/1 and p2 must be 0/0 to continue
                    if p1[0] != './.' and p2[0] != './.' and p1[0] != '0/1' and p2[0] != '0/1' and p1[0] != p2[0] and p1[0] != '0/0' and p2[0] != '1/1' and calls > minIndiv:
                        outfile.write(line)                                                         # apply count and parent filters, write out to file
            except IndexError:                                                                      # write headers to file
                if line[0] == '#':
                    outfile.write(line)
                else:
                    pass
        infile.close()
        outfile.close()
    def BestSNP(self, infile_path, outfile_path, minMinor = 1): # filter for best SNP per 150bp, minor allele count
        # written by JKK, modified by NCS
        # minMinor = minimum number of individuals carrying minor allele required to pass filter
        if os.path.isfile(self.out_dir + '/' + outfile_path) == True:                               # check if the outfile exists, give warning
            a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
            if a == 'n':
                quit()
        outfile = open(self.out_dir + '/' + outfile_path, 'w')
        infile = open(infile_path, 'r')
        scaffold = ''                                                                               # define parameters
        position = ''
        last_scaffold = ''
        last_position = 0
        cp = 0
        for line in infile:                                                                         # begin reading infile, keep best SNP of every 150 bp
            try:
                if line[0] == '#':
                    outfile.write(line)
                else:
                    scaffold = line.split('\t')[0]
                    position = int(line.split('\t')[1])                 
                    count_genos = [0, 0, 0]
                    raw_genos = [line.replace('\n','').split('\t')[i].split(':')[0] for i in range(9, self.nF2s + 11) if line.replace('\n','').split('\t')[i].split(':')[0] != './.']
                    count_genos[0] = raw_genos.count('0/0')
                    count_genos[1] = raw_genos.count('0/1')
                    count_genos[2] = raw_genos.count('1/1')
                    mincc = min(count_genos[0] + count_genos[1], count_genos[2], count_genos[1])
                    if scaffold != last_scaffold or (position - last_position) > 150:
                        if cp >= minMinor:                                                          # filter on minMinor
                            outfile.write(best_snp)
                            cp = mincc
                            last_scaffold = scaffold
                            last_position = position
                            best_snp = str(line)
                        elif mincc > cp and position != last_position:
                            cp = mincc
                            best_snp = str(line)
            except TypeError:
                pass
            except IndexError:
                pass        
        if cp >= minMinor:                                                                            # write out filtered best SNPs to file
            outfile.write(best_snp)
        infile.close()
        outfile.close()
    def ParentDepth(self, infile_path, outfile_path, p1a1_min = 100, p2a1_max = 0, p1a2_max = 0, p2a2_min = 100): # filter for allele depth of parents
        if os.path.isfile(self.out_dir + '/' + outfile_path) == True:                               # check if the outfile exists, give warning
            a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
            if a == 'n':
                quit()
        outfile = open(self.out_dir + '/' + outfile_path, 'w')
        infile = open(infile_path, 'r')
        for line in infile:
            if line[0] == '#' or line[0:3] == 'CHR':                                                # write headers to file
                outfile.write(line)
            else:                                                                                   # apply min/max filters
                p1 = line.replace('\n','').split('\t')[self.nF2s + 9].split(':')                    # amphorellae (p1)
                p2 = line.replace('\n','').split('\t')[self.nF2s + 10].split(':')                   # kunthii (p2)
                if int(p2[2].split(',')[1]) <= p2a1_max and int(p1[2].split(',')[0]) <= p1a2_max and int(p1[2].split(',')[1]) >= p1a1_min and int(p2[2].split(',')[0]) >= p2a2_min:
                        outfile.write(line)
        infile.close()
        outfile.close()
    def FilterHWE(self, infile_path, outfile_path, aa_min = 0, ab_min = 0, bb_min = 0, count = False): # filter for genotype frequencies (or counts)
        # if count == False, then aa_min, ab_min, and bb_min are all frequencies (target genotype count / called genotype count)
        # if count == True, then aa_min, ab_min, and bb_min are all integers
        infile = open(infile_path, 'r')
        if os.path.isfile(self.out_dir + '/' + outfile_path) == True:                               # check if the outfile exists, give warning
            a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
            if a == 'n':
                quit()
        outfile = open(self.out_dir + '/' + outfile_path, 'w')
        for line in infile:                                                                         # write out headers
            if line[0] == '#' or line[0:3] == 'CHR':
                outfile.write(line)
            else:                                                                                   # count genotype calls
                genotypes = [gen.split(':')[0] for gen in line.replace('\n','').split('\t')[9: self.nF2s + 9] if gen.split(':')[0] != './.']
                if count == True:
                    if genotypes.count('0/0') >= aa_min and genotypes.count('0/1') >= ab_min and genotypes.count('1/1') >= bb_min:
                        outfile.write(line)                                                         # write out the SNPs which pass the filter
                else:
                    tot = len(genotypes)
                    if genotypes.count('0/0')/tot >= aa_min and genotypes.count('0/1')/tot >= ab_min and genotypes.count('1/1')/tot >= bb_min:
                        outfile.write(line)
        infile.close()
        outfile.close()

class AdjustFiles:
    def __init__(self, out_dir = None):
        # out_dir = directory for outfiles (if None, outfiles are written to working directory)
        if out_dir != None:
            try:                                                                                    # generate directory (does not overwrite existing directory)
                self.out_dir = str(out_dir)
                os.makedirs(self.out_dir)       
            except FileExistsError:
                pass
        else:
            self.out_dir = '.'
    def MakeVCF(self, parentcalled_path, reference_path, outfile_path): # turn parentcalled VCF back into VCF with genotypes
        # parentcalled_path = file from ParentCall2
        # reference_path = any filtered VCF from pre-ParentCall2 to use as a reference for genotypes
        # outfile_path = filepath for output
        invcf = open(parentcalled_path, 'r')
        inref = open(reference_path, 'r')
        outfile = open(self.out_dir + '/' + outfile_path, 'w')
        snp_dict = {}
        for line in invcf:           
            if line[0] != '#':
                snp_dict.setdefault(line.split('\t')[0], [])
                snp_dict[line.split('\t')[0]].append(line.split('\t')[1])
        for line in inref:
            if line[0] == '#':
                outfile.write(line)
            if line.split('\t')[0] in snp_dict.keys() and line.split('\t')[1] in snp_dict[line.split('\t')[0]]:
                outfile.write(line)
        inref.close()
        invcf.close()
        outfile.close()
    def SplitVCF(self, invcf_path): # input VCF, output per-scaffold VCFs in vcfs directory
        # infile_path is the path to the main (multi-scaffold) VCF
        # generates 1 directory ('vcfs'), files assume title 'scaffold' followed by the indexed character from full chromosome title
        try:
            os.makedirs(self.out_dir + '/vcfs')
        except FileExistsError:
            pass
        infile = open(invcf_path, 'r')
        headers = []
        per_scaff = {}
        for line in infile:
            if line[0] == '#':
                headers.append(line)
            else:
                per_scaff.setdefault(line[self.scaffold_index], [])
                per_scaff[line[self.scaffold_index]].append(line)
        for i in per_scaff.keys():
            outfile = open(self.out_dir + '/vcfs/scaffold' + str(i) + '.vcf')
            for item in per_scaff[i]:
                outfile.write(item)
            outfile.close()
    def SplitFASTA(self, infasta_path): # input genome fasta, output per-scaffold fastas in fastas directory
        # infasta_path is the path to the genome file (fasta)
        # generates 1 directory ('fastas'), files assume title 'scaffold' followed by the indexed character from full chromosome title
        try:
            os.makedirs(self.out_dir + '/fastas')
        except FileExistsError:
            pass
        with open(infasta_path, 'r') as infile:
            done = []
            for line in infile:
                if line[0] == '>':
                    done.append(line[self.scaffold_index + 1])
                    if len(done) > 1:
                        outfile.close()
                    outfile = open(self.dir + '/fastas/scaffold' + str(line[self.scaffold_index]) + '.fasta', 'w')
                    outfile.write(line)                   
                else:
                    outfile.write(line)
        outfile.close()

class LepMAP3: # requires java and Lep-MAP3
    def __init__(self, out_dir = None, nF2s = 238, scaffold_index = 12, lepmap_filepath = '~/Lep-MAP3/bin', pedigree_filepath = 'pedigree.txt'):
        # out_dir = directory for outfiles (if None, outfiles are written to working directory)
        # nF2s = number of F2 individuals in the VCF
        # scaffold_index = index (base 0) within the title of each scaffold which can be used to define that scaffold
        #   this is used to name files with a single digit rather than the entire scaffold name
        #   for example, my scaffold 0 title starts with 'PGA_scaffold0...', so I only use the character in position 12 (0) to define it 
        # lepmap_filepath = filepath to Lep-MAP3 bin
        # pedigree_filepath = filepath to pedigree
        try:                                                                                        # verify that nF2s is an integer
            self.nF2s = int(nF2s)        
        except TypeError:
            print('NF2s must be an integer.')
            quit()            
        try:                                                                                        # verify that scaffold_index is an integer
            self.scaffold_index = int(scaffold_index)
        except TypeError:
            print('Scaffold index must be an integer.')
            quit()
        if out_dir != None:
            try:                                                                                    # generate directory (does not overwrite existing directory)
                self.out_dir = str(out_dir)
                os.makedirs(self.out_dir)       
            except FileExistsError:
                pass
        else:
            self.out_dir = '.'
        if os.path.isfile(lepmap_filepath) == True:                                                 # verify that input files exist
            self.lepmap_filepath = lepmap_filepath
        else:
            print('Invalid Lep-MAP3 filepath: ' + str(lepmap_filepath))
        if os.path.isfile(pedigree_filepath) == True:
            self.pedigree_filepath = pedigree_filepath
        else:
            print('Invalid pedigree filepath: ' + str(pedigree_filepath))
    def ParentCall(self, infile_path, outfile_path): # apply ParentCall2 module to VCF
        if os.path.isfile(self.out_dir + '/' + outfile_path) == True:                               # check if the outfile exists, give warning
            a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
            if a == 'n':
                quit()
        outfile = open(self.out_dir + '/' + outfile_path, 'w')
        infile = open(infile_path, 'r')
        for line in infile:                                                                         # retain headers
            if line[0] == '#':
                outfile.write(line)
        outfile.close()                                                                             # apply module
        outstring = str('java -cp ' + self.lepmap_filepath + ' ParentCall2 data=' + self.pedigree_filepath + ' vcfFile=' + infile_path + ' removeNonInformative = 1 >> ' + self.out_dir + '/' + outfile_path)
        os.system(outstring)
    def OrderMarkersGenome(self, invcf_path): # apply OrderMarkers2 module to parentcalled VCFs (per scaffold) using genome
        # generates 3 directories (vcfs, maps, orders) from input VCF (parentcalled)
        # output files assume title 'scaffold' followed by the indexed character from full chromosome title
        infile = open(invcf_path, 'r')
        try:
            for item in ['vcfs', 'maps', 'orders']:
                os.makedirs(self.out_dir + '/' + item)
        except FileExistsError:
            pass
        headers = []
        snp_dict = {}
        for line in infile:                                                                                                     # this is a parentcalled VCF
            try:
                if line[0] == '#' or line[0:3] == 'CHR':                                                                        # retain headers from VCF
                    headers.append(line)
                elif int(line[self.scaffold_index]) == int(line[self.scaffold_index]):                                          # find current scaffold
                    snp_dict.setdefault(line[self.scaffold_index], [])
                    snp_dict[line[self.scaffold_index]].append(line.split('\t')[1])                                             # add SNP position to list of scaffold's SNPs
                    if len(snp_dict[line[self.scaffold_index]]) == 1:                                                           # if this is the first SNP from a scaffold
                        if len(snp_dict.keys()) == 1:                                                                           # write out the previous one (it is finished)
                            if os.path.isfile(self.out_dir + '/vcfs/scaffold' + str(line[self.scaffold_index]) + '.vcf') == True:   
                                a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
                                if a == 'n':
                                    quit()
                            scaffold_vcf = open(self.out_dir + '/vcfs/scaffold' + str(line[self.scaffold_index]) + '.vcf', 'w')
                            for header in headers:
                                scaffold_vcf.write(header)                                   
                            if os.path.isfile(self.out_dir + '/maps/scaffold' + str(line[self.scaffold_index]) + '.txt') == True:   
                                a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
                                if a == 'n':
                                    quit()
                            scaffold_map = open(self.out_dir + '/maps/scaffold' + str(line[self.scaffold_index]) + '.txt', 'w')
                        else:
                            scaffold_vcf.close()
                            scaffold_map.close()  
                            if os.path.isfile(self.out_dir + '/vcfs/scaffold' + str(line[self.scaffold_index]) + '.vcf') == True:   
                                a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
                                if a == 'n':
                                    quit()
                            scaffold_vcf = open(self.out_dir + '/vcfs/scaffold' + str(line[self.scaffold_index]) + '.vcf', 'w')
                            for header in headers:
                                scaffold_vcf.write(header)                                
                            if os.path.isfile(self.out_dir + '/maps/scaffold' + str(line[self.scaffold_index]) + '.txt') == True:   
                                a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
                                if a == 'n':
                                    quit()
                            scaffold_map = open(self.out_dir + '/maps/scaffold' + str(line[self.scaffold_index]) + '.txt', 'w')
                        order_index = 0                                             
                    order_index += 1
                    scaffold_vcf.write(line)
                    scaffold_map.write(str(order_index) + '\n')
            except ValueError:
                print('Funky line!' + '\n' + line)
        scaffold_vcf.close()
        scaffold_map.close()
        for i in snp_dict.keys():   # apply order markers to each scaffold using genome-determined marker order
            map_path = self.out_dir + '/maps/scaffold' + str(i) + '.txt'
            vcf_path = self.out_dir + '/vcfs/scaffold' + str(i) + '.vcf'
            if os.path.isfile(self.out_dir + '/orders/scaffold' + str(i) + '.txt') == True:   
                a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
                if a == 'n':
                    quit()
            outfile = self.out_dir + '/orders/scaffold' + str(i) + '.txt'
            outstring = str('java -cp ' + self.lepmap_filepath + ' OrderMarkers2 evaluateOrder=' + map_path + ' data=' + vcf_path + ' sexAveraged=1 improveOrder=0 > ' + outfile)
            os.system(outstring) 
    def OrderMarkersLepMAP(self, invcf_path, outmap_name, iterations = 3): # apply OrderMarkers2 module to parentcalled VCF (single) without using genome
        # outmap_name = title of single map file to be generated and used for ordering 
        # iterations = number of iterations after initial ordering
        # generates 1 directory (orders), order files assume title 'scaffold' followed by the indexed character from full chromosome title
        # output order files which have been ordered multiple times (iterated) include a '.' followed by the iteration count (ex: scaffold1.3.txt is scaffold1.txt which has undergone 3 iterations of OrderMarkers2)
        infile = open(invcf_path, 'r')
        try:
            os.makedirs(self.out_dir + '/orders')
        except FileExistsError:
            pass
        map_dict = {}
        for line in infile:                                                                            # this is a parentcalled VCF
            if line[0] != '#' and line[0:3] != 'CHR':
                map_dict.setdefault(line[self.scaffold_index], 0)                                      # find current scaffold
                map_dict[line[self.scaffold_index]] += 1
        sorted_keys = sorted(map_dict.keys())
        if os.path.isfile(self.out_dir + '/' + outmap_name) == True:   
            a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
            if a == 'n':
                quit()
        map_path = self.out_dir + '/' + outmap_name       
        mapfile = open(map_path, 'w')                                                                   # generate map file for separating chromosomes
        for n, i in enumerate(sorted_keys):
            for snp in range(0, map_dict[i]):
                mapfile.write(str(n + 1) + '\n')
        mapfile.close()
        for n, i in enumerate(sorted_keys):                                                             # generate order file to be written when we call Lep-MAP3 
            outfile_path = self.out_dir + '/orders/scaffold' + str(i) + '.txt'            
            if os.path.isfile(outfile_path) == True:   
                a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
                if a == 'n':
                    quit()            
            outstring = str('java -cp ' + self.lepmap_filepath + ' OrderMarkers2 map=' + map_path + ' data=' + invcf_path + ' chromosome=' + str(n + 1) + ' sexAveraged=1 > ' + outfile_path)
            os.system(outstring)
            last_outfile = str(self.out_dir + '/orders/scaffold' + str(i))
            if iterations > 0:                                                                          # iterate over the ordering step as determined by parameter 'iterations'
                for it in range(0, iterations):
                    if os.path.isfile(self.out_dir + '/orders/scaffold' + str(i) + '.' + str(it + 1) + '.txt') == True:   
                        a = input('This file already exists! Do you want to continue? (y/n)' + '\t')
                        if a == 'n':
                            quit()
                    it_outfile = self.out_dir + '/orders/scaffold' + str(i) + '.' + str(it + 1) + '.txt'
                    outstring = str('java -cp ' + self.lepmap_filepath + ' OrderMarkers2 evaluateOrder=' + outfile_path + ' data=' + invcf_path + ' chromosome=' + str(n + 1) + ' sexAveraged=1 > ' + it_outfile)
                    os.system(outstring)
                    last_outfile = self.out_dir + '/orders/scaffold' + str(i) + '.' + str(it + 1)
    
class PlotPrep: # requires NumPy
    def __init__(self, out_dir = None, nF2s = 238, scaffold_index = 12, bin = 1000000, scaff_lengths = {'0' : 68081433, '1' : 67113551, '2' : 65712110, '3' : 57328190, '4' : 52746053, '5' : 52024118, '6' : 46780370, '7' : 43486057}):
        # out_dir = directory for outfiles (if None, outfiles are written to working directory)
        # nF2s = number of F2 individuals in the VCF
        # scaffold_index = index (base 0) within the title of each scaffold which can be used to define that scaffold
        #   this is used to name files with a single digit rather than the entire scaffold name
        #   for example, my scaffold 0 title starts with 'PGA_scaffold0...', so I only use the character in position 12 (0) to define it             
        # the scaff_lengths dictionary is specific to the P. kunthii genome assembly and should be adjusted if using a different genome (this is used for binning)     
        import numpy as np
        try:                                                                                        # verify that nF2s is an integer
            self.nF2s = int(nF2s)        
        except TypeError:
            print('NF2s must be an integer.')
            quit()            
        try:                                                                                        # verify that scaffold_index is an integer
            self.scaffold_index = int(scaffold_index)
        except TypeError:
            print('Scaffold index must be an integer.')
            quit()
        try:                                                                                        # verify that bin is an integer
            self.bin = int(bin)
        except TypeError:
            print('Bin must be an integer.')
            quit()
        if out_dir != None:
            try:                                                                                    # generate directory (does not overwrite existing directory)
                self.out_dir = str(out_dir)
                os.makedirs(self.out_dir)       
            except FileExistsError:
                pass
        else:
            self.out_dir = '.'
        self.scaff_lengths = scaff_lengths
    def GetCoords(self, invcf_path, inord_path): # input per-scaffold or single VCF directory and per-scaffold order directory, output coordinates in coordinate directory
        # invcf_path is the path to the directory of vcf files (scaffold0.vcf, scaffold1.vcf... scaffold7.vcf)
        #   these do not have to be translated back into genotypes (i.e., not the output of ParentCall2), but they DO have to be the same length as what is used in OrderMarkers2
        # inord_path is the path to the directory of order files (scaffold0.txt, scaffold1.txt... scaffold7.txt)
        # generates 1 directory ('coords'), files assume title 'scaffold' followed by the indexed character from full chromosome title
        try:
            os.makedirs(self.out_dir + '/coords')    # generate coords directory
        except FileExistsError:
            pass
        if type(invcf_path) == list:                        # define list of VCFs (can be a collection of individual VCFs or a single VCF)
            vcf_list = invcf_path
        else:
            vcf_list = [str(invcf_path)]
        if type(inord_path) == list:                        # define list of orders (output from OrderMarkers2, there should be 1 per linkage group)
            ord_list = inord_path
        else:
            print('Order files must be given in a list!')
            quit()
        bp = {}
        for vcf in vcf_list:                                # for each VCF, extract SNP positions (bp) per scaffold
            infile = open(vcf, 'r')
            for line in infile:
                if line[0] != '#' and line.split()[0] != 'CHR':
                    scaffold = line[self.scaffold_index]
                    bp.setdefault(scaffold, [])                         
                    bp[scaffold].append(line.split()[1])
            infile.close()
        cm = {}
        for ord in ord_list:                                # for each order file (output from OrderMarkers2), extract order of markers per scaffold
            scaffold = ord[8]
            cm.setdefault(scaffold, {})                
            infile = open(ord, 'r')
            for line in infile:
                if line[0] != '#' and len(line.split()) > 1:
                    cm[scaffold].setdefault(line.split()[0], line.split()[1])
            infile.close()
        for key in bp.keys():                               # make sure that each chromosome has the correct number of SNPs from each input (order and VCF)
            if len(cm[key]) != len(bp[key]):
                print('Files do not align! Too many markers!')
                quit()
            else:
                outfile = open(self.out_dir + '/coords/scaffold' + str(key) + '.txt', 'w')
                for idx, cms in cm[key].items():            # assign indexed SNP (and corresponding cM value) to basepair position, write out to file (per scaffold)
                    outfile.write(str(bp[key][int(idx) - 1]) + '\t' + str(cms) + '\n')
                outfile.close()
    def GetGCContent(self, infas_path): # input per-scaffold fasta directory, output gc-content in gc directory
        # infas_path is the path to the directory of fasta files (scaffold0.fasta, scaffold1.fasta... scaffold7.fasta)
        # generates 1 directory ('gc'), files assume title 'scaffold' followed by the indexed character from full chromosome title
        try:
            os.makedirs(self.dir + '/gc')                                                                                  # generate GC directory
        except FileExistsError:
            pass
        infastas = [file for file in os.listdir(infas_path) if 'fasta' in file.split('.') and file[0:8] == 'scaffold']       # loop through fastas, write out GC content per bin
        for fasta in infastas:
            tot = 0
            gc = 0
            n = 0
            infile = open(infas_path + '/' + fasta, 'r')
            outfile = open('gc/scaffold' + fasta[8] + '.txt', 'w')
            for line in infile:
                if line[0] != '>' and len(line.replace('\n','')) > 0:
                    for nuc in line.replace('\n',''):
                        n += 1
                        if nuc in ['G', 'C', 'g', 'c']:
                            gc += 1
                        if nuc in ['G', 'C', 'A', 'T', 'g', 'c', 'a', 't']:
                            tot += 1
                        if n % self.bin == 0:
                            outfile.write(str(round(gc/tot, 4)) + '\n')
                            tot = 0
                            gc = 0
            outfile.write(str(round(gc/tot, 4)))
            infile.close()
            outfile.close()
    def GetSNPDensity(self, invcf_path): # input per-scaffold VCF directory, output snp densities in snpden directory
        # invcf_path is the path to the directory of vcf files (scaffold0.vcf, scaffold1.vcf... scaffold7.vcf)
        #   these do not have to be translated back into genotypes (i.e., not the output of ParentCall2), but they DO have to be the same length as what is used in OrderMarkers2
        # generates 1 directory ('snpden'), files assume title 'scaffold' followed by the indexed character from full chromosome title
        invcfs = [file for file in os.listdir(invcf_path) if file[0:8] == 'scaffold' and 'vcf' in file.split('.')]
        try:
            os.makedirs(self.out_dir + '/snpden')                                                   # generate snpden directory
        except FileExistsError:
            pass
        for vcf in invcfs:                                                                          # loop through VCFs
            snp_list = []
            scaffold = vcf[8]
            infile = open(invcf_path + '/' + vcf, 'r')
            outfile = open('snpden/scaffold' + str(scaffold) + '.txt', 'w')
            for line in infile:
                if line[0] != '#' and line[0:3] != 'CHR':                                           # record scaffold lengths, SNP locations
                    scaffold_len = self.scaff_lengths[line[self.scaffold_index]]
                    snp_list.append(int(line.split('\t')[1]))
            for i in range(0, scaffold_len, self.bin):                                              # bin SNPs, write out counts
                x = [marker for marker in snp_list if marker > i and marker <= i + self.bin]
                outfile.write(str(len(x)) + '\n')
            outfile.close()
            print('Linkage Group ' + str(scaffold) + ' has ' + str(len(snp_list)) + ' SNPs')        # print per-scaffold SNP counts
    def GetRE(self, instk_path, target): # input STK and repeat element, output per-scaffold repeat element counts in directory
        # instk_path is the path to the Stockholm file generated from RepeatModeler
        # generates 1 directory (named from target RE family), files assume title 'scaffold' followed by the indexed character from full chromosome title
        try:
            os.makedirs(self.out_dir + '/' + target.lower())                                   # generate directory with name of RE family
        except FileExistsError:
            pass
        infile = open(instk_path, 'r')
        x = 0
        n_families = 0
        count = 0
        re_dict = {}    
        for line in infile:
            if line[0:7] == '#=GF TP':                                                         # find repeat element families which match target, count them
                if target in line.replace('\n','').split(';'):
                    x = 1
                    n_families += 1
                else:
                    x = 0
            if x == 1 and line[0] != '#' and len(line.split('_')) > 5:                         # extract locations, write out to file, print counts
                re_dict.setdefault(line[self.scaffold_index], [])
                re_dict[line[self.scaffold_index]].append(line.split(':')[1].split('-')[0])
                count += 1
        print(str(count) + ' ' + str(target) + ' repeat elements detected across a total of ' + str(n_families) + ' families')
        for key in sorted(re_dict.keys()):
            print('Linkage Group ' + str(key) + ' contains ' + str(len(re_dict[key])) + ' ' + str(target) + ' elements')
            with open(self.out_dir + '/' + target.lower() + '/scaffold' + str(key) + '.txt', 'w') as outfile:
                for i in range(0, self.scaff_lengths[key], self.bin):
                    x = [element for element in re_dict[key] if int(element) > i and int(element) <= i + self.bin]
                    outfile.write(str(len(x)) + '\n')
    def GetGenoFreqs(self, invcf_path, outfile_path): # input VCF (single or list), output per-scaffold and overall genotype frequencies into one output file
        # invcf_path is the path (or list of paths) to the VCF(s) to use to calculate genotype freqs. This has to have genotypes in it (cannot be raw output of ParentCall2)
        # generates single file per VCF (outfile_path) with scaffold and overall SNP counts and genotype frequencies
        if type(invcf_path) == list:
            infile_list = [str(path) for path in invcf_path]                        # arrange input into a list
        else:
            infile_list = [str(invcf_path)]
        for n, file in enumerate(infile_list):                                      # loop through each file
            infile = open(file, 'r')            
            freq_dict = {'Overall' : [0, 0, 0]}
            for line in infile:   
                if line[0] != '#' and line[0:3] != 'CHR':
                    freq_dict.setdefault(str(line[self.scaffold_index]), [0, 0, 0])
                    genotypes = [line.split('\t')[i].split(':')[0] for i in range(9, self.nF2s + 9) if line.split('\t')[i].split(':') != './.']
                    freq_dict[str(line[self.scaffold_index])][0] += genotypes.count('0/0')
                    freq_dict[str(line[self.scaffold_index])][1] += genotypes.count('0/1')
                    freq_dict[str(line[self.scaffold_index])][2] += genotypes.count('1/1')
                    freq_dict['Overall'][0] += genotypes.count('0/0')
                    freq_dict['Overall'][1] += genotypes.count('0/1')
                    freq_dict['Overall'][2] += genotypes.count('1/1')
            for i in freq_dict.keys():  # loop through each scaffold and the 'overall' data, write out to file
                outfile = open(self.out_dir + '/' + outfile_path[n])                        
                outfile.write(str(i) + '\t' + str(sum(freq_dict[i])) + '\t' + str(freq_dict[i][0]/sum(freq_dict[i])) + '\t' + str(freq_dict[i][1]/sum(freq_dict[i])) + '\t' + str(freq_dict[i][2]/sum(freq_dict[i])) + '\n')
                outfile.close()
    def PrintMetrics(self, incoord_path): # input file for one scaffold with coordinates (from 'coord' directory), print out basic metrics
        # invcf_path = filepath for large (all scaffolds) VCF
        if os.path.isfile(incoord_path) == False:                     # verify that input file exists
            print('File does not exist: ' + incoord_path)
            quit()
        x = []
        y = []
        n = 0
        g = 0
        l = 0
        with open(incoord_path, 'r') as infile:
            for line in infile:
                if line[0] != '#' and len(line.split()) > 1:
                    x.append(int(line.split()[0]))                          # x and y are used to calcualte correlation coeff
                    y.append(float(line.split()[1].replace('\n','')))
                    n += 1                                                  # n = number of SNPs
                    if float(line.split()[1].replace('\n','')) - l > g:     # g and l are used to find greatest distance between SNPs
                        g = float(line.split()[1].replace('\n','')) - l
                    l = float(line.split()[1].replace('\n',''))
        print('Correlation coeff:' + '\t' + str(np.corrcoef(x, y)[0][1]))
        print('Number of SNPs:' + '\t' + str(n))
        print('Map distance: '+ '\t' + str(max(y)) + ' cM')
        print('Greatest distance between SNPs: ' + '\t' + str(g) + ' cM')

class Plots: # requires matplotlib (pyplot, ticker, patches) and NumPy
    def __init__(self, scaff_order = None, scaff_lengths = None, reordered_lengths = None):
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        import matplotlib.patches as mpatches
        import numpy as np
        if scaff_order == None:
            self.scaff_order =  {'3' : '1', '7' : '2', '2' : '3', '4' : '4', '5' : '5', '6' : '6', '0' : '7', '1' : '8'}     # my scaffold were out of order so I fixed them using these variables
        if scaff_lengths == None:
            self.scaff_lengths =  [0, 68081433, 67113551, 65712110, 57328190, 52746053, 52024118, 46780370, 43486057]
        if reordered_lengths == None:
            self.reordered_lengths = [0, 57328190, 43486057, 65712110, 52746053, 52024118, 46780370, 68081433, 67113551]
    def Fig23(self, infile_path, outfile_path): # input frequency file (single or list), output barplot (if more than 1, arranged from top to bottom)
        # infile_path is the path (or list of paths) to the genotype frequencies to use to make barplot
        # generates barplot of genotype frequencies (outfile_path)
        if type(infile_path) == list:
            infile_list = [str(path) for path in infile_path]
        else:
            infile_list = [str(infile_path)]
        plot_list = []
        cat_plot_list = []
        for file in infile_list:
            plot_dict = {}
            infile = open(file, 'r')
            for line in infile:
                scaffold = str(line.replace('\n','').split('\t')[0])
                try:
                    plot_dict[scaffold] = [line.replace('\n','').split()[2], line.replace('\n','').split()[3], line.replace('\n','').split()[4]]
                except IndexError:
                    print('ERROR with this line: ' + line)
            plot_list.append(plot_dict)
        width = 0.2
        for n, plot_dict in enumerate(plot_list):
            plt.subplot(len(infile_list), 1, n + 1)
            aa = [float(val[0]) for key, val in plot_dict.items()]
            ab = [float(val[1]) for key, val in plot_dict.items()]
            bb = [float(val[2]) for key, val in plot_dict.items()]
            x = np.arange(len(plot_dict.keys()))
            plt.bar(x - width, aa, width, color = 'navy')
            plt.bar(x, ab, width, color = 'purple')
            plt.bar(x + width, bb, width, color = 'crimson')
            plt.ylim(0, 0.55)  
            plt.xticks(range(0, len(plot_dict.keys())), labels = ['' for tick in range(0, len(plot_dict.keys()))])
            plt.yticks([0.0, 0.125, 0.25, 0.375, 0.5], labels = ['0.0', '', '0.25', '', '0.50'])    
        plt.xticks(x, self.plot_labels)
        plt.xlabel('Linkage Groups')
        plt.savefig(self.out_dir + '/' + outfile_path)
    def Fig4(self, lepmap_coords_dir, outfile_path, reverse = ['1', '2', '3', '4', '7'], plot_this = False): # generate dotplot for comparing genome-guided and genome-independent marker ordering
        # lepmep_coords_dir = filepath for directory which includes the coordinates to use in this plot
        # reverse = scaffolds to reverse (normally not known until initial attempt)
        # plot_this = boolean, do you want to plot this graph? Otherwise, it will just calculate corrcoef.
        lepmap_files = [file for file in os.listdir(lepmap_coords_dir)]
        plot_dict = {}
        all_bps = []
        scaff_lines = []
        x = []
        y = []
        for file in lepmap_files:                                   # loop through coordinate files, extract coordinates
            old_scaffold = str(file[8])
            new_scaffold = self.scaff_order[str(file[8])]
            scaffold_length = self.scaff_lengths[int(file[8]) + 1]
            infile = open(lepmap_coords_dir + '/' + file, 'r')
            raw_coords = {}
            for line in infile:
                if line[0] != '#' and len(line.split()) > 1:
                    raw_coords.setdefault(int(line.split()[0]), float(line.split()[1].replace('\n','')))
            infile.close()        
            adjusted_coords1 = {self.scaff_lengths[int(old_scaffold) + 1] - bp : cm for bp, cm in raw_coords.items()}               # adjust raw coordinates to position relative to scaffold start position
            adjusted_coords2 = {sum(self.reordered_lengths[0:int(new_scaffold)]) + bp : cm for bp, cm in adjusted_coords1.items()}  # adjust coordinates to new position relative to beginning of genome
            adjusted_coords3 = {cm : [] for cm in adjusted_coords2.values()}                                                        # flip dictionary to order markers which share cM values within a scaffold
            for bp, cm in adjusted_coords2.items():
                adjusted_coords3[cm].append(bp)
            plot_dict.setdefault(new_scaffold, [])                                                                                  # rename scaffold with new title (new number)
            sorted_cm = sorted(adjusted_coords3.keys())                                   
            for cm in sorted_cm:
                for bp in adjusted_coords3[cm]:
                    plot_dict[new_scaffold].append(bp)
                    all_bps.append(bp)
            scaff_lines.append(sum(self.reordered_lengths[0:int(new_scaffold)]))                                                    # keep track of where dotted lines go (defines chromosome boundaries)
        s_bp = sorted(all_bps)
        n = 0
        new_scaffold_list = [int(scaff) for scaff in plot_dict.keys()]
        sorted_new_scaffs = [str(scaff) for scaff in sorted(new_scaffold_list)]
        for key in sorted_new_scaffs:
            if key not in reverse:
                for bp in plot_dict[key]:
                    y.append(bp)
                    x.append(s_bp[n])
                    n += 1
            elif key in reverse:
                for i in range(0, len(plot_dict[key])):
                    y.append(plot_dict[key][len(plot_dict[key]) - 1 - i])
                    x.append(s_bp[n])
                    n += 1
        x = range(0, 450000000, 1000000)
        y = [i for i in range(0, 172000000, 1000000)]
        for i in range(200000000, 172000000, -1000000):
            y.append(i)
        for i in range(200000000, 450000000, 1000000):
            y.append(i)
        x_arr = np.array(x)
        y_arr = np.array(y)
        print('CorrCoef: ' + str(np.corrcoef(x, y)[0][1]))
        if plot_this == True:
            plt.xlabel('Genome-Determined Position (Mbp)')
            plt.ylabel('Position After Genome-Independent Ordering (Mbp)')
            plt.scatter(x_arr, y_arr, s=0.2, c='navy')     
            for z in scaff_lines:
                plt.axhline(y = z, ls = ':', color = 'lightgrey')
                plt.axvline(x = z, ls = ':', color = 'lightgrey')
            plt.xticks(range(0, max(s_bp), 50000000), ['0', '', '10', '', '20', '', '30', '', '40', ''])
            plt.yticks(range(0, max(s_bp), 50000000), ['0', '', '10', '', '20', '', '30', '', '40', ''])
            plt.savefig(outfile_path)   
    def Fig5(self, bin, indir, include, re, outfile_path, colors = None, labels = None): # generate plot with linkage maps, snpden, gc-content, repeat element distributions
        # bin must be an integer
        # indir = directory which has the following folders to be used for this plot: coords, snpden, gc, and the one named for the target repeat element family
        # include = list of scaffolds to include in plot
        # re = name of directory within indir which includes repeat element data
        # you can adjust your own colors and labels, I included mine here
        if type(bin) != int:
            print('Bin must be an integer!')
            quit()
        if colors == None:
            colors = {'0' : 'navy', '1' : 'indigo', '2' : 'mediumvioletred', '3' : 'mediumslateblue', '4' : 'hotpink', '5' : 'deeppink', '6' : 'rebeccapurple', '7' : 'palevioletred'}
        if labels == None:
            labels = {'0' : 'LG 7', '1' : 'LG 8', '2' : 'LG 3', '3' : 'LG 1', '4' : 'LG 4', '5' : 'LG 5', '6' : 'LG 6', '7' : 'LG 2'}
        if type(include) == list:
            include_list = [str(i) for i in include]
        else:
            include_list = [str(include)]
        gc_dict = {str(scaffold) : [] for scaffold in include_list}                                           # get GC content per scaffold, save to dictionary
        gc_paths = [filepath for filepath in os.listdir(indir + '/gc') if str(filepath[8]) in include_list]   # this will become a heatmap
        for gc_file in gc_paths:
            infile = open(indir + '/gc/' + gc_file, 'r')
            for line in infile:
                if len(line) > 0:
                    gc_dict[str(gc_file[8])].append(float(line.replace('\n','')))
            infile.close()
        re_dict = {str(scaffold) : [] for scaffold in include_list}                                             # get repeat element content per scaffold, save to dictionary
        re_paths = [filepath for filepath in os.listdir(indir + '/' + re.lower()) if str(filepath[8]) in include_list]   # this will become a heatmap
        for re_file in re_paths:
            infile = open(indir + '/' + re.lower() + '/' + re_file, 'r')
            for line in infile:
                if len(line) > 0:
                    re_dict[str(re_file[8])].append(float(line.replace('\n','')))
            infile.close()
        sd_dict = {str(scaffold) : [] for scaffold in include_list}                                             # get SNP density per scaffold, save to dictionary
        sd_paths = [filepath for filepath in os.listdir(indir + '/snpden') if str(filepath[8]) in include_list]   # this will become a histogram
        for sd_file in sd_paths:
            infile = open(indir + '/snpden/' + sd_file, 'r')
            for line in infile:
                if len(line) > 0:
                    sd_dict[str(sd_file[8])].append(float(line.replace('\n','')))
            infile.close()
        cm_dict = {str(scaffold) : [] for scaffold in include_list}                                             # get cM values per scaffold, save to dictionary
        bp_dict = {str(scaffold) : [] for scaffold in include_list}                                             # get basepair positions per scaffold, save to dictionary
        co_paths = [filepath for filepath in os.listdir(indir + '/coords') if str(filepath[8]) in include_list]   # together, these become the scatterplots
        for co_file in co_paths:
            infile = open(indir + '/coords/' + co_file, 'r')
            for line in infile:
                if len(line) > 0:
                    bp_dict[str(co_file[8])].append(float(int(line.replace('\n','').split()[0])/bin))
                    cm_dict[str(co_file[8])].append(float(line.replace('\n','').split()[1]))
            infile.close()
        fig, ax = plt.subplots(2 + (2 * len(include_list)), sharex = True, gridspec_kw = {'height_ratios': [1, 0.3] + ([0.1, 0.1] * len(include_list))})
        xlim = float(max([self.scaff_lengths[str(scaffold)] for scaffold in include_list])) / bin                 # define x limit as the maximum of the lengths of scaffolds (over bin)
        patches = []
        steps = 10
        for n, scaff in enumerate(include_list):
            ax[0].scatter(bp_dict[str(scaff)], cm_dict[str(scaff)], c = self.cl[str(scaff)], alpha = 0.2, s = 1)  # generate plots
            ax[0].yaxis.set_major_locator(ticker.MultipleLocator(50))
            ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(25))
            ax[1].hist([range(0, len(sd_dict[str(scaff)]), 1)], weights = sd_dict[str(scaff)], histtype = 'step', color = self.colors[str(scaff)], alpha = 0.7)
            ax[1].xaxis.set_minor_locator(ticker.MultipleLocator(5))
            ax[1].set_xticks(*[range(0, (8 * steps), steps)], ['0', '10', '20', '30', '40', '50', '60', '70'])
            ax[1].set_xlim(-0.5, xlim + 1)
            ax[(2*n) + 2].imshow([gc_dict[str(scaff)]], cmap = 'plasma', aspect = 'auto', vmin = 0.3000, vmax = 0.4000)
            ax[(2*n) + 3].imshow([re_dict[str(scaff)]], cmap = 'viridis', aspect = 'auto', vmin = 0, vmax = 100)
            ax[(2*n) + 2].set_yticks([])
            ax[(2*n) + 3].set_yticks([])
            patches.append(mpatches.Patch(color = self.colors[str(scaff)], label = self.labels[str(scaff)]))
        ax[0].legend(handles = patches)
        plt.savefig(outfile_path)
    def Fig5_Scatterplots(self, indir, outfile_path, labels = None, scaff_order = None):
            # indir = directory which has the following folders to be used for this plot: coords, snpden, gc, and the one named for the target repeat element family
            # include = list of scaffolds to include in plot
            # re = name of directory within indir which includes repeat element data
            # you can adjust your own colors and labels, I included mine here
            if labels == None:
                labels = {'0' : 'LG 7', '1' : 'LG 8', '2' : 'LG 3', '3' : 'LG 1', '4' : 'LG 4', '5' : 'LG 5', '6' : 'LG 6', '7' : 'LG 2'}
            if scaff_order == None:
                scaff_order = ['3', '7', '2', '4', '5', '6', '0', '1']
            co_paths = [filepath for filepath in os.listdir(indir + '/coords')]   # together, these become the scatterplots
            cm_dict = {}
            bp_dict = {}
            for co_file in co_paths:
                cm_dict.setdefault(str(co_file[8]), [])                                         # get cM values per scaffold, save to dictionary
                bp_dict.setdefault(str(co_file[8]), [])                                         # get basepair positions per scaffold, save to dictionary
                infile = open(indir + '/coords/' + co_file, 'r')
                for line in infile:
                    if len(line) > 0:
                        bp_dict[str(co_file[8])].append(int(line.replace('\n','').split()[0]))
                        cm_dict[str(co_file[8])].append(float(line.replace('\n','').split()[1]))
                infile.close()
            print(cm_dict)
            print(bp_dict)
            fig, ax = plt.subplots(4, 2, sharex = True, sharey = True)
            patches = []
            n = 0
            for scaff in scaff_order:
                if n % 2 == 0:
                    ax[int(n/2), 0].scatter(bp_dict[str(scaff)], cm_dict[str(scaff)], c = 'darkblue', alpha = 0.2, s = 2)  # generate plots
                    ax[int(n/2), 0].set_xticks([0, 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000], ['0', '10', '20', '30', '40', '50', '60', '70'])                    
                    n += 1
                else:
                    ax[int((n - 1)/2), 1].scatter(bp_dict[str(scaff)], cm_dict[str(scaff)], c = 'darkblue', alpha = 0.2, s = 2)  # generate plots
                    n += 1
            ax[0, 0].yaxis.set_major_locator(ticker.MultipleLocator(50))
            ax[0, 0].yaxis.set_minor_locator(ticker.MultipleLocator(25))       
            ax[0, 1].yaxis.set_major_locator(ticker.MultipleLocator(50))
            ax[0, 1].yaxis.set_minor_locator(ticker.MultipleLocator(25))              
            plt.savefig(outfile_path)


