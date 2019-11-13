#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:39:11 2019

@author: hannes
"""



class data:

    def __init__(self, name):
            self.name = name
            self.df = {}
            self.settings = {}


#initiate class
data = data('Data') 

#import os
import sys
import re
import getopt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
from multiprocessing import Pool,cpu_count



#instructions
form='\nanotomat.py\n\n\n\t\
                --genome <genome.fasta>\n\t\
                --pos <start_positions.txt> *exact format see below*\n\t\
                --cov <coverage> *file with \'#chr bp coverage\' from "samtools depth -a -Q 0 <file.bam>"*\n\t\
                --mincov <int> minimal coverage of the START pos for a reannotation*default=0*\n\t\
                --mincov_exon <int> minimal coverage for annotation of exon after intron*default=0*\n\t\
                --gt use less stringency in defining introns *G._.G defines intron*(default=GT_AG)\n\t\
                --cores <int> number of cores used *default=(cpu_count-1)*\n\t\
                --name <gene name> *default=\'genes\'*\n\t\
                --out <name for output files> *default=\'genes\'*\n\n\n\
EXAMPLE for start_positions.txt:\n\n\
Bgt_chr-01	49207	~BgtAcSP-31373	#make a comment\n\
Bgt_chr-01	284818	~BgtE-20066	#tab separated???\n\
Bgt_chr-01	300256	~BgtE-20114	#special sign \"~ \"before gene-name!!!\n\
                '



def print_anot_txt():
    txt = open('{}_anot.txt'.format(data.settings['out']), 'w+')
    for i in data.df['anot']:
        temp = ' '.join(i)
#            print('_'.join(i[:2]), 'similar', data.df['pd_met'].check.values[0])
        if '_'.join(i[:2]) in data.df['pd_met'].check.values.tolist():
#                print('anot debugging:', '_'.join(i), 'in pd_met')
            temp += ' ' + '~' + data.df['pd_met'].loc[data.df['pd_met'].check == '_'.join(i[:2]), 'gene'].values[0]
            if  len(data.df['pd_met'].loc[data.df['pd_met'].check == '_'.join(i[:2]), 'family'].values.tolist()) > 0:
                temp += ' ' + '#' + data.df['pd_met'].loc[data.df['pd_met'].check == '_'.join(i[:2]), 'family'].values[0]                
        print(temp, file=txt)
    return

def print_errors():    
    #print ERRORS
    df = data.df['faulty'].loc[data.df['faulty'].no_start == 'yes'].copy()
    if len(df) > 0:
        print('\n\n!!! Faulty genes without start codon: !!!\n\n\n')
        for i in range(len(df)):
            gene = df.iloc[i,0]
            print(gene)

    df = data.df['faulty'].loc[data.df['faulty'].no_stop == 'yes'].copy()
    if len(df) > 0:
        print('\n\n!!! Faulty genes without stop codon: !!!\n\n\n')
        for i in range(len(df)):
            gene = df.iloc[i,0]
            print(gene)
    return


def fasta_print(dic):
    txt = open('{}.fasta'.format(data.settings['out']), 'w+')
    #original seq?
    for i in dic:
        seq = dic[i]
        #header
        print(i)
        print(i, file=txt)
        #print seq
        for j in range(len(seq)//50):
            print(seq[j*50:(j+1)*50:1])
            print(seq[j*50:(j+1)*50:1], file=txt)
        if len(seq)%50 != 0:
            print(seq[(len(seq)//50)*50::])
            print(seq[(len(seq)//50)*50::], file=txt)
        print()
    return

def fasta_aa_print(dic):
    txt = open('{}_aa.fasta'.format(data.settings['out']), 'w+')
    #original seq?
    for i in dic:
        seq = dic[i]
        #header
        print(i)
        print(i, file=txt)
        #print seq
        for j in range(len(seq)//50):
            print(seq[j*50:(j+1)*50:1])
            print(seq[j*50:(j+1)*50:1], file=txt)
        if len(seq)%50 != 0:
            print(seq[(len(seq)//50)*50::])
            print(seq[(len(seq)//50)*50::], file=txt)
        print()
    return


def gff_print(dic):
    txt = open('{}.gff'.format(data.settings['out']), 'w+')
    print('##gff-version 3')
    print('##gff-version 3', file=txt)
    #original seq?
    for i in dic:
        
#            seq = dic[i][0]
        if len(dic[i][1]) > 0:
            print('###')
            print('###', file=txt)
            

# =============================================================================
#                 # print gene
# =============================================================================
            temp=[]
            lyst = dic[i][0].copy()
            lyst[2] = 'gene'
            if lyst[6] == '+':
                lyst[3] = dic[i][1][0][0]
                lyst[4] = dic[i][1][-1][1]
            else:
                lyst[3] = dic[i][1][-1][0]
                lyst[4] = dic[i][1][0][1]
            temp.append(lyst)

            temp[0][5] = '.'
            temp[0][7] = '.'
            #add exon number
            temp[0][-1] = temp[0][-1].split('.mrna1.exon')[0]  + ';Name' + \
                          temp[0][-1].split('Name')[1].split(';Parent')[0] + ';DE=NULL'
            print('\t'.join(map(str, temp[0])))
            print('\t'.join(map(str, temp[0])), file=txt)


# =============================================================================
#                 # print mrna
# =============================================================================
            temp=[]
            lyst = dic[i][0].copy()
            lyst[2] = 'mRNA'
            if lyst[6] == '+':
                lyst[3] = dic[i][1][0][0]
                lyst[4] = dic[i][1][-1][1]
            else:
                lyst[3] = dic[i][1][-1][0]
                lyst[4] = dic[i][1][0][1]
            temp.append(lyst)

            temp[0][5] = '.'
            temp[0][7] = '.'
            #add exon number
            temp[0][-1] = temp[0][-1].split('.exon')[0]  + ';Name' + \
                          temp[0][-1].split('Name')[1].split('.mrna1;Target')[0]
            print('\t'.join(map(str, temp[0])))
            print('\t'.join(map(str, temp[0])), file=txt)

# =============================================================================
#                 # print exon
# =============================================================================
            temp=[]
            for j in range(len(dic[i][1])):
                lyst = dic[i][0].copy()
                lyst[2] = 'exon'
                lyst[3] = dic[i][1][j][0]
                lyst[4] = dic[i][1][j][1]
                temp.append(lyst)
            cnt_first=0
            cnt_second=0
            for j in range(len(temp)):
                temp[j][7] = '.'
                cnt_first = cnt_second +1
                cnt_second = cnt_first + (int(temp[j][4]) - int(temp[j][3]))
                temp[j][-1] = temp[j][-1].split(' ')[0] + ' ' + str(cnt_first) + ' ' + str(cnt_second) + ' ' + '+'
                #add exon number
                temp[j][-1] = temp[j][-1].split('exon')[0] + 'exon' + str(j+1) + ';Name' + \
                              temp[j][-1].split('Name')[1]
                print('\t'.join(map(str, temp[j])))
                print('\t'.join(map(str, temp[j])), file=txt)


# =============================================================================
#                 # print CDS
# =============================================================================

            temp=[]
            left=0
            for j in range(len(dic[i][1])):
                lyst = dic[i][0].copy()
                lyst[3] = dic[i][1][j][0]
                lyst[4] = dic[i][1][j][1]
                lyst[7] = str((3-left)%3)#str(left)
                left = (left + (int(lyst[4]) - int(lyst[3]) +1 ))%3
                temp.append(lyst)
            cnt_first=0
            cnt_second=0
            for j in range(len(temp)):
                cnt_first = cnt_second +1
                cnt_second = cnt_first + (int(temp[j][4]) - int(temp[j][3]))
                temp[j][-1] = temp[j][-1].split(' ')[0] + ' ' + str(cnt_first) + ' ' + str(cnt_second) + ' ' + '+'
                #add exon number
                temp[j][-1] = temp[j][-1].split('exon')[0] + 'cds' + str(j+1) + ';Name' + \
                              temp[j][-1].split('Name')[1]
                temp[j][2] = 'CDS'
                print('\t'.join(map(str, temp[j])))
                print('\t'.join(map(str, temp[j])), file=txt)
            print('')
            print('', file=txt)
    return


def assemble_start_codon(chromo, pos):
            first_int = pos

            if data.df['genome'][chromo][int(first_int)-1:first_int+2] == 'ATG':
                print('orientation: +')
                data.settings['orient'] = '+'
            elif data.df['genome'][chromo][int(first_int)-3:first_int] == 'CAT':
                print('orientation: -')
                data.settings['orient'] = '-'
            else:
                
                print('ERROR: no start CODON!', chromo, pos)
                gene = data.df['pd_met'].loc[data.df['pd_met'].check == chromo + '_' +str(pos), 'gene'].values[0]
                data.df['faulty'].loc[data.df['faulty'].gene == gene, 'no_start'] ='yes'
            return

def assemble_seq():
    dic = {}
    dic_gff = {}
    cnt=0
    first = 0
    header=''
    seq=''
    for lyst in data.df['anot']:
            gene = data.df['pd_met'].loc[data.df['pd_met'].check == lyst[0] + '_' + lyst[1], 'gene'].values[0]
            if gene != '':
                name = gene
                
            else:
                name = '{}{:04d}'.format(data.settings['name'], cnt)
                data.df['pd_met'].loc[data.df['pd_met'].check == lyst[0] + '_' + lyst[1], 'gene'] = name
            cnt+=1
            if first == 0:
                first = 1
            else:
                dic[header] = seq
                seq=''
            header = '>'+ name + ' ' + lyst[0] + ' ' + lyst[1] + ' ' + lyst[-1]
            dic_gff[name] = [[lyst[0], data.settings['genome'], \
                    'CDS', 'start', 'stop', '100' , \
                    '+', '.', 'ID={}.mrna1.exon;Name={};Parent={}.mrna1;Target={}'.format(name,name, name, name)]]
            dic_gff[name].append([])
            lyst_no = lyst[1:]
            #sort to ascending order
            lyst_no = sorted(map(int, lyst_no))
            #splice sequences
            for i in range(len(lyst_no)):
                if i%2 == 0:
                    chromo = lyst[0]
                    first_int = lyst_no[i]
                    sec_int = lyst_no[i+1]  
                    seq += data.df['genome'][chromo][int(first_int)-1:int(sec_int)]
                    #backwards settings, store CDS sequences in descending order
                    if int(lyst[-1]) > int(lyst[1]):
                        dic_gff[name][1].append([first_int, sec_int])
                    else:
                        first_int = lyst_no[-i-2]
                        sec_int = lyst_no[-i-1]
                        dic_gff[name][1].append([first_int, sec_int])
                        dic_gff[name][0][6] = '-'
    #last contig
    if seq != '':
        dic[header] = seq
    seq=''
    header = ''

    #complement-reverse backward strand sequences
    for i in dic:
        seq = Seq(dic[i], generic_dna)
        if i.split()[2] > i.split()[3]:
            dic[i] = str(seq.reverse_complement())
    #create protein dic
    dic_aa={}
    for i in dic:
        seq = Seq(dic[i], generic_dna)
        dic_aa[i] = seq.translate()
    return dic, dic_aa, dic_gff



def annotator_multi(core):
    temp=[]
    chromo_in_use = ''
    for chromo in data.settings['chromo_lyst'][core]:
        print()
        print('core:\t', core, sep='')
        print('chromo:\t', chromo, sep='')
        print()
        if chromo in data.df['met']:
            for pos in data.df['met'][chromo]:
                print('\nannotating:', chromo, pos)
                assemble_start_codon(chromo, pos)
                
                gene = data.df['pd_met'].loc[data.df['pd_met'].check == chromo + '_' +str(pos), 'gene'].values[0]
                #skip genes without start
                if data.df['faulty'].loc[data.df['faulty'].gene==gene, 'no_start'].values[0] == 'yes':
                    continue
                #create subset coverage for chromosome
                if chromo_in_use != chromo:
                    df = data.df['cov'].loc[data.df['cov']['chr']==chromo].reset_index(drop=True)
                    coverages=df.iloc[:,1].values.tolist()
                chromo_in_use = chromo
                #skip genes without expression
                if data.settings['mincov'] != '':
                    if coverages[pos-1] < data.settings['mincov']:#df.iloc[pos-1,1] < data.settings['mincov']:
                        print('skiped lowcov:\t', chromo, pos)
                        continue
    #                    print(df)
    #                    print('small df')
    
                lyst = [chromo, str(pos)]
    
    # =============================================================================
    #                 fORWARDS STRAND!!!
    # =============================================================================
                if data.settings['orient'] == '+':
                    seq = 'ATG'
    #                    print(data.df['genome'][chromo][pos-1:pos+20])
                    cnt = pos+1
                    cov =df.iloc[cnt,1]
                    #avoid division by zero
                    if cov==0:
                        cov=10e-6
                    state='exon'
                    while True:
                        cnt+=1
                        #break if leaving chromosome
                        if cnt == len(df)-3:
    
                            data.df['faulty'].loc[data.df['faulty'].gene == gene, 'no_stop'] = 'yes'
                            break
                        new_cov =df.iloc[cnt,1]
                        #avoid division by zero
                        if new_cov==0:
                            new_cov=10e-6
    #                        print('\npos:', cnt+1) #debug
    #                        print('new_cov:', new_cov) #debug
    #                        print('cov:', cov) #debug
    #                        if state == 'exon' and (new_cov/cov) < data.settings['ratio']:# and new_cov <= data.settings['ratio_to_start']*df.iloc[pos,1]:
                        if state == 'exon' and (new_cov/cov) < data.settings['ratio'] and cov > data.settings['ratio_to_start']*df.iloc[pos,1]:
# =============================================================================
#                             if data.df['genome'][chromo][cnt:cnt+2] == 'GT' :
#     #                                print('INTRON', chromo, pos) #debug
#                                 lyst.append(str(cnt))
#                                 state = 'intron'
# =============================================================================
                            if data.settings['gt'] != '':
                                if data.df['genome'][chromo][cnt:cnt+1] == 'G':
                                    lyst.append(str(cnt))
                                    state = 'intron'
                            else:
                                 if data.df['genome'][chromo][cnt:cnt+2] == 'GT':
                                    lyst.append(str(cnt))
                                    state = 'intron'
    #                            else: #debug
    #                                print('intron candidate', chromo, cnt) #debug
    #                                print(data.df['genome'][chromo][cnt:cnt+2]) #debug
                        elif state == 'intron' and cov/new_cov <= (data.settings['ratio']) and new_cov >= data.settings['ratio_to_start']*df.iloc[pos,1] and new_cov >= data.settings['mincov_exon']:
    #                            print('EXON', chromo, pos) #debug
                            if data.settings['gt'] != '':
                                if data.df['genome'][chromo][cnt-1:cnt] == 'G':
                                    lyst.append(str(cnt+1))
                                    state = 'exon'
                            else:
                                 if data.df['genome'][chromo][cnt-2:cnt] == 'AG':
                                    lyst.append(str(cnt+1))
                                    state = 'exon'
                        if state == 'exon':
                            seq += data.df['genome'][chromo][cnt]
                        cov = new_cov
                        if len(seq)%3 == 0:
    #                            print(seq[-3:]) #debug
                            if  seq[-3:] == 'TAG' or\
                                seq[-3:] == 'TAA' or\
                                seq[-3:] == 'TGA':
                                lyst.append(str(cnt+1))
                                temp.append(lyst)
    #                                print('>'+chromo, pos, data.settings['orient']) #debug
    #                                print(seq) #debug
                                break
    
    # =============================================================================
    #                 BACKWARDS STRAND!!!
    # =============================================================================
                elif data.settings['orient'] == '-':
                    seq = 'TAC'
    #                    print(data.df['genome'][chromo][pos-21:pos:][::-1])
                    cnt = pos-3
                    cov =df.iloc[cnt,1]
                    #avoid division by zero
                    if cov==0:
                        cov=10e-6
                    state='exon'
                    while True:
                        cnt-=1
                        if cnt == 1:
                            data.df['faulty'].loc[data.df['faulty'].gene == gene, 'no_stop'] = 'yes'
                            break
                        new_cov =df.iloc[cnt,1]
                        #avoid division by zero
                        if new_cov==0:
                            new_cov=10e-6
    #                        print('\npos:', cnt+1)
    #                        print('new_cov', new_cov)
    #                        print('cov', cov)
    #                        print(data.df['genome'][chromo][cnt-1:cnt+1])
    #                        if state == 'exon' and new_cov/cov < data.settings['ratio']:# and new_cov <= data.settings['ratio_to_start']*df.iloc[pos,1]:
                        if state == 'exon' and (new_cov/cov) < data.settings['ratio'] and cov > data.settings['ratio_to_start']*df.iloc[pos,1]:
                            if data.df['genome'][chromo][cnt-1:cnt+1] == 'AC' :
    #                                print('INTRON', chromo, cnt+2)
                                lyst.append(str(cnt+2))
                                state = 'intron'
    #                            else:
    #                                print('intron candidate', chromo, cnt+2)
    #                                print(data.df['genome'][chromo][cnt+1:cnt+3])
                        elif state == 'intron' and cov/new_cov <= (data.settings['ratio']) and new_cov >= data.settings['ratio_to_start']*df.iloc[pos,1]:
                            if data.df['genome'][chromo][cnt+1:cnt+3] == 'CT':
    #                                print('EXON', chromo, cnt+1)
                                lyst.append(str(cnt+1))
                                state = 'exon'
    #                            else:
    #                                print('Exon candidate', chromo, cnt+1)
    #                                print(data.df['genome'][chromo][cnt+1:cnt+3])
    #                                sys.exit()
                        if state == 'exon':
                            seq += data.df['genome'][chromo][cnt]
                        cov = new_cov
    
                        if len(seq)%3 == 0:
    #                            print(seq[-3:])
                            if  seq[-3:] == 'ATC' or\
                                seq[-3:] == 'ATT' or\
                                seq[-3:] == 'ACT':
    #                                print(seq[-3:])
                                lyst.append(str(cnt+1))
                                temp.append(lyst)
    #                                print('>'+chromo, pos, data.settings['orient'])
                                seq = Seq(seq, generic_dna)
                                seq = str(seq.complement())
    #                                print(seq)
                                break

#                            
    return temp

def prepare_multiprocess():
    print('\n***preparing for multi-processing***\n')
    data.settings['chromo_lyst'] = {}
    for i in range(data.settings['cores']):
        data.settings['chromo_lyst'][i] = []
    lyst = []
    for chromo in data.df['genome']:
        lyst.append(chromo)
    #sort for better distribuition
    lyst = sorted(lyst, key=lambda x: len(data.df['genome'][x]), reverse=True)
    cnt = 0    
    for chromo in lyst:
        data.settings['chromo_lyst'][(cnt)%data.settings['cores']].append(chromo)
        cnt+=1
    return

def sort_dic(dic_gff, dic_aa):
    print('\n***sorting final output after multi-processing***\n')
    gff_lyst = []
    aa_lyst = []
    dic_aa_old = dic_aa
    for i in dic_gff:
#            print(i)
        temp = dic_gff[i]
        temp.append(i)
        gff_lyst.append(temp)
        temp = ''
    for i in dic_aa:
        print(i)
        print(dic_aa[i])
#            temp = dic_aa[i]
        temp = [i.split()]#, temp]
        aa_lyst.append(temp)
        temp = ''
    dic_gff = {}
    dic_aa = {}
    print(gff_lyst)
    print(aa_lyst[0])
    gff_lyst = sorted(gff_lyst, key=lambda x: (x[0][0], int(x[1][0][0])))
    aa_lyst = sorted(aa_lyst, key=lambda x: (x[0][1], int(x[0][2])))
    for i in gff_lyst:
        print('gff:', i)
        dic_gff[i[-1]] = i[:-1]
    for i in aa_lyst:
        print('aa:', ' '.join(i[0]))
        dic_aa[' '.join(i[0])] = dic_aa_old[' '.join(i[0])]
    return dic_gff, dic_aa






def main(argv):
    #set max coverage ratio for intron/exon change
    data.settings['ratio'] = .5
    #set minimal ratio for exon to initial coverage
    data.settings['ratio_to_start'] = .1
    #set variables
    data.settings['gt'] = ''
    data.settings['met'] = ''
    data.settings['name'] = 'genes'
    data.settings['genome'] =  ''
    data.settings['out'] = 'genes'
    data.settings['mincov'] = ''
    data.settings['mincov_exon'] = 0
    data.settings['cores'] = int(cpu_count()-1)
    data.df['anot'] = []
    data.df['genes'] = {}
    data.df['faulty'] ={}

    try:
       opts, args = getopt.getopt(argv,"h",["gt", "genome=", "pos=", "cov=", "mincov=", "mincov_exon=", "cores=", "name=", "out="])
    except getopt.GetoptError:
       print ('{}'.format(form))
       sys.exit()
    for opt, arg in opts:
       if opt == '-h' or opt == '-help' or opt == '--help':
          print ('{}'.format(form))
          sys.exit()
       elif opt == '--gt':
          data.settings['gt'] = 1
       elif opt == '--mincov_exon':
          data.settings['mincov_exon'] = int(arg)
       elif opt == '--genome':
          data.settings['genome'] = arg
       elif opt == '--pos':
          data.settings['met'] = arg
       elif opt == '--cov':
          data.settings['cov'] = arg
       elif opt == '--mincov':
          data.settings['mincov'] = int(arg)
       elif opt == '--cores':
          data.settings['cores'] = int(arg)
       elif opt == '--name':
          data.settings['name'] = arg
       elif opt == '--out':
           data.settings['out'] = arg
    if argv == []:
          print ('{}'.format(form))
          sys.exit()

    print()
    print(argv)
    print()






    def read_fasta():
        print('\n*** reading genome to memory ***\n')
        t=0
        seq=''
        dic={}
        ass=open('{}'.format(data.settings['genome']))
        for line in ass:
            line = line.strip('\n')
        
            #check if .fasta header
            match= re.search("^(>)", line)
        
            #first contig  
            if match and t == 0:
                header = line.split()[0].replace('>', '')
                t=1
            #if not first contig
            elif match and t==1:
        
                #add last contig, start new one
                dic[header] = seq
                header = line.split()[0].replace('>', '')
                seq=''
            else:
                seq += line
        
        #last contig
        dic[header] = seq
        seq=''
        t=0
        ass.close()
        #make all uppercase
        for i in dic:
            dic[i] = dic[i].upper()
        data.df['genome'] = dic
        print('\n*** finished reading genome ***\n')
        return 


    def read_met():
        dic = {}
        txt = open(data.settings['met'], 'r+')
        for line in txt:
            if line.strip() != '' and line[0] != '#':
                line = line.strip().split()
                if line[0] not in dic:
                    dic[line[0]]=[]
                dic[line[0]].append(int(line[1]))
        data.df['pd_met'] = pd.read_csv(data.settings['met'], '\t', header=None)
        data.df['pd_met'].columns = ['chr', 'pos', 'gene', 'family']
        data.df['pd_met']['gene'] = data.df['pd_met']['gene'].str[1:]
        data.df['pd_met']['family'] = data.df['pd_met']['family'].str[1:]
        data.df['pd_met']['check'] = data.df['pd_met'].chr.astype(str) + '_' + data.df['pd_met'].pos.astype(str)
        del data.df['pd_met']['chr']
        del data.df['pd_met']['pos']
        data.df['met'] = dic
        #create error df
        data.df['faulty'] = data.df['pd_met'].copy()
        data.df['faulty']['no_start'] = 'no'
        data.df['faulty']['no_stop'] = 'no'
        data.df['faulty'] = data.df['faulty'].loc[:,['gene', 'no_start', 'no_stop']]
        return


    def read_cov():
        print('\n*** reading coverage to memory ***\n')
        df = pd.read_csv(data.settings['cov'], '\t', header=None)
        df.columns=['chr', 'bp', 'cov']
        df['bp'],df['cov'] = df['bp'].astype(int),df['cov'].astype(int)
        del df['bp']
        print(df)
        data.df['cov']=df
        print('\n*** finished reading coverage to memory ***\n')
        return



# =============================================================================
#     !!! execution of functions !!!
# =============================================================================

    read_met()
    read_fasta()
    read_cov()

#main condition necessary for multi-threading
if __name__ == "__main__":
    main(sys.argv[1:])

    prepare_multiprocess()
    p = Pool(data.settings['cores'])
    temp = p.map(annotator_multi, list(range(data.settings['cores'])))
    for i in temp:
        data.df['anot'].extend(i)
    data.df['anot'] = sorted(data.df['anot'], key=lambda x: (x[0], int(x[1])))
    dic, dic_aa, dic_gff = assemble_seq()
    dic_gff, dic_aa = sort_dic(dic_gff, dic_aa)

    #print information in terminal
    print('\n\n\nDNA.fasta\n\n\n')
    fasta_print(dic)
    print('\n\n\nAA.fasta\n\n\n')
    fasta_aa_print(dic_aa)
    print('\n\n\nanot.gff\n\n\n')
    gff_print(dic_gff)
    print_anot_txt()
    print_errors()


            
    print('\n\nOUTFILES:\n\n\n{}.gff\n\n{}_aa.fasta\n\n{}.fasta\n\n{}_anot.txt\n'.format(data.settings['out'],\
                                                                      data.settings['out'],\
                                                                      data.settings['out'],\
                                                                      data.settings['out']))
    sys.exit()