#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 20:42:14 2022
@author: michaelgaunt
"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from rnalaura import Rnatricks as Rnat
from multiprocessing import Pool, cpu_count
from collections import deque


class Biotricks():
    def __init__(self,shuffledict={}, mc=0, seqdict={},ab=[None,None]):
        self.seqdict = seqdict
        self.ab = ab
        self.shuffledict = shuffledict
        self.mc = mc
        self.align = MultipleSeqAlignment(deque([]))

    def reformShuffle(self):
        tmp = {}
        processes_count = (cpu_count()-1)
        pool = Pool(processes_count)
        poolout = pool.map_async(rShufdict, [[shkey,self.shuffledict] for shkey in self.shuffledict.keys()],chunksize=3)
        iniRnat = Rnat('','',pool)
        pooloutGet = iniRnat.pooladmin(poolout)
        for stuff in pooloutGet:
            tmp.update(stuff)
        self.seqdict = tmp
        if self.mc:
            return self
        return tmp

    def bioobj (self, chunk = 5):
        processes_count = (cpu_count()-1)
        pool2 = Pool(processes_count)
        poolout2 = pool2.map_async(parallelBioSeq, [[ke,va] for ke, va in self.seqdict.items()], chunksize = chunk)
        iniRna = Rnat('','',pool2)
        pooloutGet2 = iniRna.pooladmin(poolout2)
        return pooloutGet2

    def bioloop(func):
        def wrapper(self):
            rnaSeqobj = []
            for ke,va in self.seqdict.items():
                rnaSeqobj.append(func(self,ke,va))
            return self.makealign(rnaSeqobj)
        return wrapper

    def makealign(self, rnaSeqobj):
        for record in rnaSeqobj:
            record.id = str(record.description)
            self.align.append(record)
        return self.align

    @bioloop
    def bioHdf5obj(self, ke, va):
        return SeqRecord(Seq(str(va.decode('utf-8')))[self.ab[0]:self.ab[1]], id= str(ke), name='', description=str(ke))

    @bioloop
    def bioobj2(self,ke, va):
        return SeqRecord(Seq(va[self.ab[0]:self.ab[1]]), id= str(ke), name='', description=str(ke))
    
    def __getitem__(self, item):
        return dict.__getitem__(self, item) % self
    
def parallelBioSeq (items):
    seqID = items[0]
    seqSeq = items[1]
    record = SeqRecord(Seq(seqSeq), id= str(seqID), name='', description=str(seqID))
    record.id = str(record.description)
    return record

def rShufdict (hand):
        collectIt = {}
        for card in hand[1][hand[0]]:
            collectIt.update({card[0]:card[1]})
        return collectIt
     