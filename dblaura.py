#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 19:10:43 2022

@author: michaelgaunt
"""
from qclaura import Qctricks as Qct
from rnalaura import Rnatricks as Rnat
import gc, re, io, json
from pathlib import Path

result_list = []

def resetit():
    global result_list
    result_list = []
    
class Dbtricks():
    def __init__(self, pathfile):
        self.pathfile=pathfile

    def subit(self, replace):
        return re.sub('(\.[^.]+?$){1,5}',replace,str(self.pathfile))

class Pangodict():

    def pangodictmake(self, resultList):
        z = re.compile(r'_')
        pangodicti = {}
        for record in resultList:
            zList = z.split(record)
            ncbiKey = zList.pop(0)
            ptype = ''.join(zList[0])
            if re.search('None',ptype):
                ptype = re.sub('None','NA\.0',ptype)
            if re.search(r'[A-Z]+\.[0-9]+\.[0-9]+\.[0-9]+$', ptype):
                ptypeTruncate = re.sub(r'([A-Z]+\.[0-9]+\.?[0-9]*)\.?[0-9]*$', r'\1', ptype)
            else:
                ptypeTruncate = ptype
            pangodicti[ncbiKey] = {'pangotype': [ptype, ptypeTruncate]}
            if len(zList) == 1:
                continue
            else:
                country = ''.join(zList[1])
                ptime = zList[2]
                if len(zList) > 3:
                    species = zList[3]
                    pangodicti[ncbiKey].update({'species': species})
                pangodicti[ncbiKey].update({'country': country})
                pangodicti[ncbiKey].update({'ptime': ptime})
        return pangodicti

class ParellelFAio(Dbtricks, Pangodict):
    def __init__(self, pathfile, pangotime, nreps = 10, qc = 0, write = 0, rangeloc =list([0,550]),chunk=2):        
        self.pathfile = pathfile
        self.pangotime = pangotime
        self.nreps = nreps
        self.qc = qc
        self.write = write
        self.rangeloc = rangeloc
        self.pangolisti = []
        super().__init__(pathfile)

    def _grep(self, pattern, rootfile):
        c = re.compile(pattern)
        with io.open(rootfile, "r", encoding="utf-8") as fin:
            if c.search(fin.read()):
                return True
            else:
                return False

    def _choices(self):
        try:
            if str(self.qc) == '1':
                fout = Path(self.subit('.qc.fa'))
            elif int(self.pangotime) == 0 and str(self.nreps) == '0':
                fout = Path(self.subit('.0.fa'))
            elif int(self.pangotime) == 1 or int(self.pangotime) == 2:
                fout = Path(self.subit('.v' + str(self.pangotime) + '.pango.fa'))
        except Exception as error:
            print ('The QC or Pangotype flags are not aligned {}'.format(error) )
            return False
        return fout

    def parallDictWrit(self, delIDs=[]):
        if self.pangotime == 0 and self.nreps > 0 and self.write == 0:
            return {}
        if self._grep('^>([^_]+?_){2,8}.*', self.pathfile) == False and self._grep(' ', self.pathfile) == False and self._grep('__', self.pathfile) == False:
            print (f"Please check {self.pathfile} it may be already parsed")
            return 
        fout = self._choices()
        if fout.is_file():
            fout.unlink()
        iniQct = Qct(self.pathfile)
        x, y, z, w2, y2, z2, z3, z4 = self._compiles()
        w = iniQct.wcompile()
        # print (f'Look here for DB {fout}')
        pool = iniQct.makepool()
        for record in iniQct.getalign():
            pool.apply_async(parallelwrite, args=(record,delIDs,fout,w,w2,x,y,y2,z,z2,z3,z4, self.pangotime, self.write, self.qc, self.nreps,  self.rangeloc),  callback = self.log_result)
        inirna = Rnat('','',pool)
        inirna.pooladmin()
        del pool, iniQct, inirna
        gc.collect()
        # jsonout = re.sub (r'fa$','v2.json', str(fout))
        if str(self.nreps) == '0':
            resetit()
            return fout
        if str(self.pangotime) == '1' or str(self.pangotime) == '2':
            if not str(self.nreps) == '0':
                pangodicti = super().pangodictmake(result_list)
                # with open(Path(jsonout), "w") as djout:
                #     djout.write(str(pangodicti))
                return pangodicti
            else:
                print(f"Replication is {self.nreps} so no genetic type or epi calculation permitted for pangotime {self.pangotime}")

    def log_result(self,result):
        result_list.append(result)

    def _compiles(self):
        x = re.compile(r':')
        y = re.compile(r'[ ]')
        z = re.compile(r'\|')
        w2 = re.compile(r'(^[^\.^_]+?)\.?[0-9]*?(_.*)')
        y2 = re.compile(r'[_]{2,4}')
        z2 = re.compile(r'^([^_]+?_[^_]+?)_.*')
        # TODO: Change to 5 for species
        z3 = re.compile(r'^(([^_]+?_){4}).*')
        z4 = re.compile(r'_$')
        return x, y, z, w2, y2, z2, z3, z4

def parallelwrite(record, delIDs, fout, w, w2, x, y, y2, z, z2, z3, z4, pangotime, write, qc, nreps, rangeloc):
    if len(delIDs) > 0 and str(qc) == '1':
        for recDel in delIDs:
            tmpP = w.sub(r'\1', recDel)
            if re.search(tmpP, str(record.id)):
                return ''
    headerFa = x.sub('', str(record.description))
    headerFa = y.sub('_', headerFa)
    headerFa = z.sub(':', headerFa)
    headerFa = y2.sub('_', headerFa)
    headerFa = w2.sub(r'\1\2', headerFa)
    if str(qc) == '1':
        headerFa = re.sub(r'$',r'%%', headerFa)
    elif str(nreps) == '0':
        headerFa = w.sub(r'\1', headerFa)
    elif str(pangotime) == '1':
        headerFa = z2.sub(r'\1', headerFa)
    elif str(pangotime) == '2':
        headerFa = z3.sub(r'\1', headerFa)
        headerFa = z4.sub('', headerFa)
    if str(qc) == '1' or str(nreps) == '0' or str(write) == '1':
        with open (fout, 'a') as fw:
            fw.write('>%s\n%s\n' % (headerFa, str(record.seq)[rangeloc[0]:rangeloc[1]]))
    if str(pangotime) == '1' or str(pangotime) == '2':
        return headerFa
    else:
        return False

'''
if __name__ == '__main__':
    startfile = 'SARScov2_5UTR_1-0.25mil.aln.fa'
    pangotime = 2
    reps = 1
    qc = 0
    pathway = '/Volumes/data/5UTR_millcuarto/testTrue1'
    startfile = 'SARScov2_5UTR_1-0.25mil.aln.fa'
    startpath = Path(pathway, startfile)
    file = Path(pathway, startfile)
    delIDs = []
    if re.search(r'[12]', str(pangotime)) or str(qc) == '1' or str(reps) == '0':
        iniParse = ParellelFAio(startpath, pangotime, reps, qc)
        pangodicti = iniParse.parallDictWrit(delIDs)
    print (pangodicti)
'''