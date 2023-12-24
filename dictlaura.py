#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:46:20 2022

@author: michaelgaunt
"""
from __future__ import print_function
import random

def sortDict (mystuff):
    return sorted(mystuff.items(), key=lambda x: x[1]['freq'], reverse=True)
    
def shuffleDict (d):
    skeys = list(d.keys())
    random.shuffle(skeys)
    return {k:d[k] for k in skeys}

if __name__ == '__main__':

    mystuff = {'OL369093': {'pangotypeTrun': ['B.1.1', 'B.1.1'], 'pangotype': ['B.1.1.7', 'B.1.1.7'], 'freq': 2, 'country': ['USA', 'USA'], 'ptime': ['2021-04-29', '2021-04-29']}, 'OL968587': {'pangotypeTrun': ['AY.100'], 'pangotype': ['AY.100'], 'freq': 1, 'country': ['USA'], 'ptime': ['2021-12-01']}, 'OL902478': {'pangotypeTrun': ['NA\\.0'], 'pangotype': ['NA\\.0'], 'freq': 1, 'country': ['USA'], 'ptime': ['2021-09-14']}, 'MT831889': {'pangotypeTrun': ['A.1'], 'pangotype': ['A.1'], 'freq': 1, 'country': ['USA'], 'ptime': ['2020-06-16']}, 'OL902472': {'pangotypeTrun': ['NA\\.0'], 'pangotype': ['NA\\.0'], 'freq': 4, 'country': ['USA'], 'ptime': ['2021-09-14']}, 'MT831889': {'pangotypeTrun': ['A.1'], 'pangotype': ['A.1'], 'freq': 1, 'country': ['USA'], 'ptime': ['2020-06-16']}}
    sortedD = sortDict(mystuff)
    shuffleDs = shuffleDict(mystuff)
    print (sortedD)
