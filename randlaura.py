#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 14:55:28 2022

@author: michaelgaunt
"""
import random
import numpy as np
from pathlib import Path
import pandas as pd


def dot_bracket_to_pairs(ss_string):
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
    rangelist = list(range(1,len(sorted_df)))
    isorted_df = sorted_df.reset_index(drop=True)
    return isorted_df

def makestems53(df):
    last5,last3 = 0, 0
    indexNew = 0
    stems5  = dict()
    stems3 = dict()
    stems5[indexNew] = []
    stems3[indexNew] = []
    flag = 0
    for index, col in df.iterrows():
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
    stemDict5 = _refineStems(stems5)
    stemDict3 = _refineStems(stems3)
    return stemDict5, stemDict3

def _refineStems (stems):
    stemDictTmp = {k: v for k, v in stems.items() if len(v) > 0}
    stemDict = {i+1: v for i, v in enumerate(stemDictTmp.values())}
    return stemDict

def shuffle_parts(array, slices):
    count = 0
    aout = np.array([])
    for s in slices:
        if count == 1:
            array = np.flipud(a)
            random.Random(6).shuffle(array[slice(*s)])
            aout = np.flipud(array)
        else:
            random.Random(6).shuffle(a[slice(*s)])
            count += 1
    print (aout)
    
if __name__ == '__main__':

    size = 1750
    chunk = 1
    startfile = 'SARScov2_5UTR_1-0.25mil.aln.fa'
    # startfile = 'SARScov2_5UTR_1-0.25mil.reduced.fa'
    pathway = '/Volumes/data/5UTR_millcuarto/testTrue2'
    rnaref = 'rnaref.fa'
    tracker = Path(pathway, '.laura')
    locusVals = ''
    # a = np.arange(100)
    # info = Counter()
    db_string = '((((((((((..(((((((.......)))))))......).((((((.......))))))..)))))))))'
    # sequence =  list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrs')
    sequence = list('012345678-----------------------------------------------------876543210')
    a = np.asarray(sequence)
    df = dot_bracket_to_pairs(db_string)
    stemDict5, stemDict3 = makestems53(df)
    # print(stemDict3)
    start = stemDict5[1][0]
    end = start + int(len(stemDict5[1])-1)
    start2 = stemDict3[1][-1]
    end2 = start2 + int(len(stemDict3[1])-1)
    print (end2)
    shuffle_parts(array=a, slices=((start,end), (start, end)))
        
    # print(stemDict3)
    
    # newSeq = []
    # stemcount = 0
    # for x in range(0,len(sequence)):
    #     if x not in stemDict:
    #         newSeq.append(x)
    #     else:
    #         newSeq.append(stemDict[stemcount])
    #         stemcount += 1
    #         continue
