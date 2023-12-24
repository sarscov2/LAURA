#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:26:22 2022
This combines 
* ParallelHDF5 = converting to HDF5 format
* ParellelFAio = writing and header formating
* RnaQC = removing noisy sequences
@author: michaelgaunt
"""

from rnalaura import Rnatricks as Rnat
from biolaura2 import Biotricks as Bt
from multiprocessing import Pool, cpu_count
import re, gc
from pathlib import Path
from Bio import SeqIO
from Bio import AlignIO
import h5py
from collections import Counter, defaultdict


class Qctricks():
    def __init__(self, pathfile = '', ab=[None,None], chunk=1, seqdict={}):
        self.pathfile=pathfile
        self.seqdict = seqdict
        self.ab = ab

    def makepool(self):
        processes_count = (cpu_count()-1)
        return Pool(processes_count)
        
    def wcompile (self):
        return re.compile (r'^([^\.^_]+?)[._].*')

# TODO: this needs removing and redirecting to dblaura
    def subit(self, replace):
        return re.sub('(\..*){1,3}',replace,str(Path(self.pathfile).name))

    def getalign(self):
        return AlignIO.read(Path(self.pathfile),'fasta')

class ParallelHDF5(Qctricks):
    def __init__(self, pathfile, ab=[None,None]):
        # self.pathfile = pathfile
        self.fout = pathfile
        # self.fout = Path (Path(self.pathfile).parent, self.subit('.hdf5'))
        self.seqdict = {}
        self.ab = ab
        super().__init__(pathfile, ab, 1, {})

    def paraFAtoHDF5(self):
        w = self.wcompile()
        pool = self.makepool()
        if self.fout.is_file():
            self.fout.unlink()
        for record in self.getalign():
            pool.apply_async(parallelHDF5, args=(record, self.fout, w))
        admin = Rnat('','',pool)
        admin.pooladmin()

    def paraHDF5(self, align):
        w = self.wcompile()
        # pool = self.makepool()
        if self.fout.is_file():
            self.fout.unlink()
        hd5path = '5UTR/'
        for record in align:
            header = w.sub(r'\1',str(record.description))
            with h5py.File(self.fout, 'a') as h5file:
                h5file[hd5path + header] = str(record.seq)

    def openHDF5(self):
        if self.fout.is_file():
            myhdf5 = h5py.File(self.fout, 'r')
            item = myhdf5['5UTR']
            # print(list(myhdf5.keys()))
            keys = item.keys()
            for key in keys:
                self.seqdict[key] = item[key][()]
            iniBt = Bt({},1,self.seqdict, self.ab) 
            return iniBt.bioHdf5obj()

class ParellelFAio(Qctricks):
    def __init__(self, pathfile, rangeloc =list([0,550]),chunk=2, fastaAlign={}):        
        self.rangeloc = rangeloc
        self.pathfile = pathfile
        self.fastaAlign = fastaAlign
        super().__init__(pathfile, rangeloc)
        
    def _removeDels(self, delIDs, rangeloc):
        for record in self.fastaAlign:
            if str(record.id) in delIDs:
                print ("Removed {}".format(record.id))
            else:
                record.seq = record.seq[rangeloc[0]-1: rangeloc[1]]
                yield record

    def parallRepwrit(self, delIDs=[], rangeloc = [None,None]):
        fout = re.sub(r'fa$', 'qc.fa', str(self.pathfile))
        with open (Path(fout), 'w') as f:
            SeqIO.write(self._removeDels(delIDs, rangeloc), f, "fasta")
        return fout
                
class RnaQC(ParellelFAio, Qctricks):
    def __init__(self, filepath, rangeloc = [None,None], chunk=2, fastaAlign = {}):
        # self.fullgenome = AlignIO.read(Path(path, fullgenome),'fasta')
        self.filepath = filepath
        self.fastaAlign = fastaAlign
        self.chunk = chunk
        self.ambigDict = defaultdict(int)
        super().__init__(filepath, rangeloc, chunk, fastaAlign)

    def _ambigCounts(self,poolout, x):
        genomeFreqs = {ids:freqs for ncbi in poolout.get() for ids, freqs in ncbi.items()}
        for nucid, freqs in genomeFreqs.items():
            for k,v in freqs.items():
                if not x.match(k):
                    self.ambigDict[nucid] += v
                else:
                    self.ambigDict[nucid] += 0

    def coreQC (self, threshold = 25, locusVals = '1-400'):
        rangeloc = [int(i) for i in re.split('-', locusVals)]
        self.fastaAlign = self.getalign()
        x = re.compile (r'[ATGC-]')
        pool = self.makepool()
        poolout = pool.map_async(
            rnaQCcount, [record for record in  self.getalign()],chunksize=self.chunk)   
        admin = Rnat('','',pool)
        admin.pooladmin()
        self._ambigCounts(poolout, x)
        del pool
        ambigSorted = dict(sorted(self.ambigDict.items(), key=lambda item: item[1]))
        delIDs = [ids for ids, ambig in ambigSorted.items() if ambig > threshold] 
        ambigList = [aSv for aSv in ambigSorted.values()]
        del self.ambigDict
        fout = self.parallRepwrit(delIDs,rangeloc)
        gc.collect()
        return fout, delIDs, ambigSorted, ambigList
    
def parallelwrite(record,delIDs,fout,w,w2,x,y,y2,z,rangeloc):
    header = w.sub(r'\1',str(record.description))
    flag = 0
    for di in delIDs:
        delid = w.sub(r'\1',di)
        if header == delid:
            flag = 1
        else:
            continue
    if flag == 1:
        return 0
    headerFa = x.sub('', str(record.description))
    headerFa = y.sub('_', headerFa)
    headerFa = z.sub(':', headerFa)
    headerFa = w2.sub(r'\1\2', headerFa)
    headerFa = y2.sub(r'_', headerFa)
    headerFa = re.sub(r'$',r'%%', headerFa)
    with open (fout, 'a') as fw:
        fw.write('>%s\n%s\n' % (headerFa, str(record.seq)[rangeloc[0]:rangeloc[1]]))

def parallelHDF5(record,fout,w):
    hd5path = '5UTR/'
    header = w.sub(r'\1',str(record.description))
    with h5py.File(fout, 'a') as h5file:
        # h5file[hd5path + header] = str(record.seq)
        h5file[hd5path + header] = str(record.seq)

def rnaQCcount(record):
    seqFreqs = {}
    listN = list(str(record.seq).upper())
    tmp = {k:v for k,v in Counter(listN).items()}
    seqFreqs.update({str(record.id):tmp})
    return seqFreqs

if __name__ == '__main__':
    path = '/Volumes/Data/05UTR_millcuarto/testTrue1/phylip314'
    # path = '/Volumes/Data/5UTR_millcuarto/larder'
    genFile = 'SARScov2_5UTR_1-0.25mil.aln.315.1.hdf5'
    # alignFile = 'SARScov2_5UTR_1-0.75mil.aln.fa '
    testout = 'testout3.hdf5'
    fullgenome = Path(path, genFile)
    ref = 'rnaref.fa'
    rnaref = Path(path, ref)
    fullgenome = Path(path, genFile)
    iniHdf5 = ParallelHDF5(fullgenome)
    stuff = iniHdf5.openHDF5()
    fout = Path(path,testout)
    iniHdf5w =  ParallelHDF5(fout)
    iniHdf5w.paraHDF5(stuff)
    
    # iniHdf5 = ParallelHDF5(fout)
    # print (iniHdf5.openHDF5())
    
    # inipara = RnaQC(fullgenome)
    # delIDs, ambigSorted, ambigList = inipara.countNs()
    # print (len(ambigList))
