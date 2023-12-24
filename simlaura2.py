#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 14:55:28 2022

@author: michaelgaunt
"""

from random import Random, randint, seed, random
import numpy as np
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from rnalaura2 import Rnatricks as Rnat
import re, gc

class Simwrite():
    def __init__(self, record, pathway, myseed, fileout, simout, stemList, newname):
        self.record = record
        self.pathway = pathway
        self.myseed = myseed
        self.fileout = fileout
        self.simout = simout
        # This is stemDict5 and stemDict3
        self.stemList = stemList
        self.newname = newname
    
    def writesim (self, reps):
            loopNo = [k for k in self.stemList[0].keys()]
            with open (self.simout, 'w') as fout:
                for loopN in loopNo * reps: 
                    randomNo = randint(0,10000)
                    fout.write(
                        '>{}_loop{}_random={}_seed={}\n{}\n'.format
                        (self.newname, loopN, randomNo, self.myseed, self.loopShuffle 
                         (randomNo, loopN))
                        )
    
    def writestems(self):
        combinedDict = list(self.common_entries(self.stemList[0], self.stemList[1]))
        with open (Path(
                self.pathway, 'stem_positions_' + self.newname + '.txt'), 'w') as fout2:
            fout2.write('Stem{0}5\' end{0}3\' end\n'.format("\t"))
            for cd in combinedDict:
                trueLoc5 = [c+1 for c in cd[1]]
                trueLoc3 = [c+1 for c in cd[2]]
                fout2.write(
                    "{}\t{}-{}\t{}-{}\n".format(
                    cd[0],min(trueLoc5),max(trueLoc5),min(trueLoc3),max(trueLoc3)
                    ))
                
    # https://stackoverflow.com/questions/16458340/python-equivalent-of-zip-for-dictionaries
    def common_entries(self, *dcts):
        if not dcts:
            return
        for i in set(dcts[0]).intersection(*dcts[1:]):
            yield (i,) + tuple(d[i] for d in dcts)
            
    def loopShuffle(self, randNo, loopN):
        a = np.asarray(list(self.record.seq))
        start = self.stemList[0][loopN][0]
        end = start + int(len(self.stemList[0][loopN])-1)
        # end = start + int(len(stemDict5[loopN]))
        start2 = self.stemList[1][loopN][-1]
        end2 = start2 + int(len(self.stemList[1][loopN]))
        revStem5 = len(self.record.seq) - end2
        revStem3 = revStem5 + (end2 - start2 - 1)
        # revStem3 = revStem5 + (end2 - start2)
        return ''.join(self._shuffle_parts(a, randNo, slices=((start,end), (revStem5,revStem3))))
    
    def _shuffle_parts(self, a, randNo, slices):
         count = 0
         aout = np.array([])
         for s in slices:
             if count == 1:
                 arraytmp = np.flipud(a)
                 Random(randNo).shuffle(arraytmp[slice(*s)])
                 aout = np.flipud(arraytmp)
             else:
                 Random(randNo).shuffle(a[slice(*s)])
                 count += 1
         return aout    

    def __repr__(self):
         return repr(self.stemList)
     
    def __str__(self):
        return str(self.stemList)

class Simclean(Simwrite):
    def __init__(self, record, pathway, myseed, fileout, simout, isorted_df, stemList = '', newname = ''):
        self.record = record
        self.pathway = pathway
        self.myseed = myseed
        self.fileout = fileout
        self.simout = simout
        self.isorted_df = isorted_df
        self.stemList = stemList
        self.newname = newname
        super().__init__(record, pathway, myseed, fileout, simout, stemList, newname)

    def _makestems53(self):
        last5,last3, indexNew = 0, 0, 0
        stems5, stems3 = dict(), dict()
        stems5[indexNew], stems3[indexNew] = [], []
        flag = 0
        for index, col in self.isorted_df.iterrows():
            if col["5\'"] - last5 == 1 and last3 - col["3\'"] == 1:
                # print (col["5\'"])
                if flag == 0:
                    stems5[indexNew].append(col["5\'"] - 1)
                    stems3[indexNew].append(col["3\'"] + 1)
                stems5[indexNew].append(col["5\'"])
                stems3[indexNew].append(col["3\'"])
                flag = 1
            else:
                indexNew +=1
                stems5[indexNew] = []
                stems3[indexNew] = []
                flag = 0
            last5 = col["5\'"]
            last3 = col["3\'"]
        self.stemDict5 = self._refineStems(stems5)
        self.stemDict3 = self._refineStems(stems3)
        del self.isorted_df
        return self
    
    def _refineStems (self, stems):
        stemDictTmp = {k: v for k, v in stems.items() if len(v) > 0}
        stemDict = {i+1: v for i, v in enumerate(stemDictTmp.values())}
        return stemDict
        
    def _loopClean1(self):
        stemList = list()
        for stemDict in [self.stemDict5, self.stemDict3]:
             tmp = {k:v for k,v in stemDict.items() if any(i < 400 for i in v )}     
             stemList.append(tmp)
        return stemList[0], stemList[1]
    
    def _loopClean2(self,stemDict1, stemDict2):
        return {k:stemDict2[k] for k in stemDict1.keys() if k in stemDict2}
    
    def _loopClean3(self, stemDict5ii,stemDict3ii):
        self.stemList = list()
        for sd in [stemDict5ii, stemDict3ii]:
            tmp = {k:v for k,v in sd.items() if max(v) - min(v) > 3}
            self.stemList.append(tmp)
        return self
        
    def _loopCleanK (self):
        stemListtmp = list()
        for sd in [self.stemList[0], self.stemList[1]]:
            tmpii = {i+1:v for i,v in enumerate(sd.values())}
            stemListtmp.append(tmpii)
        del self.stemList
        self.stemList = stemListtmp
        return self
    
    def stemDriver(self):
        self._makestems53()
        stemDict5i, stemDict3i = self._loopClean1()
        stemDict3ii = self._loopClean2(stemDict5i, stemDict3i)
        stemDict5ii = self._loopClean2(stemDict3ii, stemDict5i)
        self._loopClean3(stemDict5ii,stemDict3ii)
        return self._loopCleanK()
    
class Shufflesim (Simclean):
    def __init__ (self,record='', pathway='', myseed='', simout='', fileout='', isorted_df=''):
        self.record = record
        self.pathway = pathway
        self.myseed = myseed
        self.simout = simout
        self.fileout = fileout
        self.isorted_df = isorted_df
        super().__init__(record, pathway, myseed, fileout, simout, isorted_df)
    
    def rnastrucout (self):
        self.newname = re.sub(r"(^[^_]+)_.*",r"\1", str(self.record.id))
        fileout = Path(self.pathway,'rnastruc_' + self.newname + '.dot')
        if not fileout.is_file():
            with open (fileout, 'w') as fout:
                fout.write('>{}\n{}\n{}\n'.format(str(self.record.id), str(self.record.seq), rnastruct))
        return self

    def dot_bracket_to_pairs(self, ss_string):
        '''
        Takes dot and bracket string, returns dataframe
        with paired bases.
        If any invalid characters are in the structure, it
        will interpret them as dots, as it only reads parentheses.
        '''
        index_list = []
        pairs = {}
        for index, char in enumerate(ss_string):
            if char == '(':
                index_list.append(index)
            if char == ')':
                try:
                    # pair to last item in the list in dictionary
                    pairs[index_list.pop()] = index
                except IndexError:
                    print(f'Invalid structure, found extra \')\' in position {index}')
    
        if len(index_list) != 0:
            for item in index_list:
                print(f'Invalid structure, found extra in \'(\' in position {item}')
    
        df_pairs = pd.DataFrame(pairs.items(), columns=["5\'", "3\'"])
        sorted_df = df_pairs.sort_values(by=['5\''], ascending=True)
        # rangelist = list(range(1,len(sorted_df)))
        self.isorted_df = sorted_df.reset_index(drop=True)
        return self
    
    # def __str__(self):
    #     return str(self.isorted_df)
    
    # def __repr__(self):
    #     return repr(self.isorted_df)
    
    @staticmethod
    def getsimout(pathway, rnaref):
        records = SeqIO.parse(Path(pathway,rnaref),'fasta')
        simname = '_'.join([re.sub(r"(^[^_]+)_.*",r"\1", r.id) for r in records])
        simout = Path(pathway, 'stemshuffle.sim.' + simname + '.fa')
        if simout.is_file():
            simout.unlink()
        return simout
            
if __name__ == '__main__':

    pathway = '/Volumes/data/5UTR_millcuarto/testTrue2'
    rnaref = 'rnaref.fa'    
    reps = 1
    myseed = 10
    records = SeqIO.parse(Path(pathway,rnaref),'fasta')

    stem1 = {'1':[[1,2,3,4,5],[14,13,12,11,10]]}
    x = [[x,y]for x,y in zip(stem1['1'][0],stem1['1'][1])]
    randloc = round(random()*len(stem1['1'][0]))-1
    x[randloc].reverse()
    stem2 = {1:[z for z in zip(*x)]}
    print (stem2)
    
    raise SystemExit()
    iniSs = Shufflesim()
    simout = iniSs.getsimout(pathway, rnaref)
    del iniSs
    rnaname, rnastruct,rnaseq = '','',''
    for record in records:
        iniRn = Rnat(record.seq, pathway)
        rnastruct, score = iniRn.getStruct()
        seed(myseed)
        iniSs = Shufflesim(record, pathway, myseed, simout)
        # Need the method below to generate the outfile name for
        # positioning and simulation
        iniSs.rnastrucout()
        iniSs.dot_bracket_to_pairs(rnastruct)
        iniSs.stemDriver()
        iniSs.writesim()
        iniSs.writestems()
        myseed  += 1 
    