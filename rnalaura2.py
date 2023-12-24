#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 22 22:07:17 2022
There are 3 separate algorithms here:
    * Rnaloop
    * Alignstruc
    * Plotting
@author: Michael Gaunt
"""

from Bio import AlignIO, SeqIO
from pathlib import Path
from multiprocessing import Pool, cpu_count
import re, sys, shutil
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import numpy as np
from itertools import islice
from subprocess import Popen, PIPE
import deepdish as dd

result_list = []

def resetit():
    global result_list
    result_list = []

class Findrnafold():
    def __init__(self, rnaref='', pathDir='', pool=''):
        self.pathDir = pathDir
        self.rnaref = rnaref
        self.rnaStruc = ''
        self.pool = pool
        
    def tutranl (self):
        TUtran = {84:85}
        return str(self.rnaref).translate(TUtran)

    def rnafoldTest(self, locate):
        for w in locate:
            if Path(w).is_file():
                rnafold = Path(w)
                break
        check = Popen(
            ['echo "AA" | %s' % (rnafold)], stdout=PIPE, 
            stderr=PIPE, shell=True)
        stdout, stderr = check.communicate()
        if stdout.decode().split('\n')[0] == 'AA':
            return rnafold
        else:
            return False

    def errorout(self):
        return "Please insert the full path to 'rnafold' under the '-rp' flag.\n Alternatively, place it in either '~/bin', '/usr/local/bin'.\n"
    
    def errorout2(self):
        return "At for OSX at 'https://www.tbi.univie.ac.at/RNA/download/osx/macosx/ViennaRNA-2.5.0-MacOSX.dmg'\n or for Ubuntu ''https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_22_04/viennarna_2.5.1-1_amd64.deb' \n or for other Linux at 'https://www.tbi.univie.ac.at/RNA/#download'\n."

    def rnaprog (self, user='~'):
        locate = shutil.which('rnafold')
        locate2 = shutil.which('RNAfold')
        if not Path(str(locate)).is_file() or not Path(str(locate2)).is_file():
            for loc in ['~/bin', '/usr/local/bin', user]:
                if not locate2:
                    break
                binpath = Path(loc).resolve()
                for prog in ['RNAfold', 'rnafold']:
                    locate = Path(binpath, prog)
                    try:
                        if locate.is_files():
                            locate2 = False
                            break
                        else:
                            raise Exception (
                                "I don't think 'rnafold' is installed\n." 
                                + "If it is installed then:" 
                                + self.errorout
                                + "'rnafold can be installed with:\n"
                                + "'conda install -c bioconda viennarna', or\n"
                                + self.errorout2()
                                + "Finally, there is a copy of 'rnafold' via,\n"
                                + "docker run -it --rm /phylogenetics/laura" 
                                +" --entrypoint /bin/bash"
                                )
                            sys.exit(0)
                    except FileNotFoundError as fnf_error:
                        print(fnf_error)
            rnafold = self.rnafoldTest([locate])
        else:            
            rnafold = self.rnafoldTest([locate, locate2])
        if rnafold:
            return rnafold
        else:
            raise Exception (
                "Can find 'rnafold' at {}\n".format(rnafold)
            + "but it just ain't working for me.\n"
            + self.errorout()
            )
            sys.exit(0)
                
class Rnatricks (Findrnafold):
    def __init__(self, rnaref='', pathDir='', pool=''):
        self.pathDir = pathDir
        self.rnaref = rnaref
        self.rnaStruc = ''
        self.pool = pool

    def tutranl (self):
        TUtran = {84:85}
        return str(self.rnaref).translate(TUtran)

    def getStruct(self,rnafold='/usr/local/bin/RNAfold'):
        # rnafold = self.rnaprog()
        # /usr/local/bin/RNAfold
        process = Popen(['echo %s | %s' % (self.tutranl(),rnafold)],
                        stdout=PIPE, 
                        stderr=PIPE,
                        shell=True
                        )
        stdout, stderr = process.communicate()
        rnaStrucTmp = stdout.decode().split('\n')
        rnaStruc = rnaStrucTmp[1].split(' ')
        score = float(rnaStruc[1].lstrip('(').rstrip(')'))
        return rnaStruc[0], score
    
    def writeIt(self):
        rnaStrucRef = Path(self.pathDir,'rnaref2.dot')
        if rnaStrucRef.is_file():
            rnaStrucRef.unlink()
        with open (Path(self.pathDir, 'rnaref2.dot'), 'w') as fout:
            fout.write("%s\n%s" % (str(self.rnaref), self.rnaStruc))

    def pooladmin(self, poolout=''):
        self.pool.close()
        self.pool.join()
        self.pool.terminate()   
        if poolout:
            return poolout.get()
        
    def testit (self):
        return "this worked"

class Rnaloop (Rnatricks):
    def __init__(self, pathDir, rnaref):
        super().__init__(rnaref, pathDir)
        # This maybe better expressed using AOP
        # print (Path(pathDir, rnaref))
        self.rnaref = re.sub(r'-', '', str(SeqIO.read(Path(pathDir, rnaref), 'fasta').seq))
        self.pathDir = pathDir
        self.score = 0
        self.rnaStruc = ''

    def _deloopM (self, z, a, g, g2):
        groupsF = z.search(self.rnaStruc)
        groupsR = a.search(self.rnaStruc)
        try:
            countF = groupsF.group(g).count('(')
        except:
            countF = 0
        try:
            countR = groupsR.group(g2).count(')')
        except:
            countR = -1
        return countF, countR

    def deloopMega(self):
        self.rnaStruc, self.score = self.getStruct()
        # super().writeIt()
        z = re.compile(r'^([\.]{0,50})([\(]{2,8})([\.]{2})')
        z2 = re.compile(r'^([\.]{0,50})([\(]{2,8}[\.]{0,2}[\(]{2,8}[\.]{0,10}[\)]{2,8}\.*[\)]{2,8})([\.]{1,8})([\(]{2,8})([\.]{1,8})')
        a = re.compile(r'([\.]{1})([\)]{1,8}\.?[\)]{1,8})([\.]{0,50})$')
        countF, countR = self._deloopM (z, a, 2, 2)
        if countF == countR:
            insertF = '?' * countF
            self.rnaStruc = z.sub(r'\1%s\3'%(insertF),self.rnaStruc)
            self.rnaStruc = a.sub(r'\1%s\3'%(insertF),self.rnaStruc)
        else:
            countF2, countR2 = self._deloopM (z2, a, 4, 2)
            if countF2 == countR2:
                insertF2 = '?' * countF2
                self.rnaStruc = z2.sub(r'\1\2\3%s\5'%(insertF2),self.rnaStruc)
                self.rnaStruc = a.sub(r'\1%s\3'%(insertF2),self.rnaStruc)

    def _clearStems (self, x, x2, term):
        xvar = x.findall(str(self.rnaStruc))
        try:
            for var in xvar:
                insertion = term * var.count('.')
                # rnaStruc = re.sub(x2 %(var), r'\1%s\2'%(insertion), rnaStruc)
                self.rnaStruc = re.sub(x2 %(var), r'\1%s?\2'%(insertion), self.rnaStruc)
        except:
            pass

    def _gapLoc(self):
        return [m.start() for m in re.finditer('-', self.rnaref)]

    def gapcorrect (self, loopLocs):
        result = list()
        gaps = self._gapLoc()
        for rangeLoop in loopLocs:
            for x in rangeLoop:
                for y in gaps:
                    if y <= x:
                        x = x + 1
                result.append(x)
        return  [(result[i],result[i+1]) for i in range(0,len(result),2)]

    def _loopInsert(self):
        x = re.compile(r'[\)]{4}([\.]{1,5})[\(]{1}[\.]{1}[\(]{4}')
        x2 = r'([\)]{4})%s[\(]{1}([\.]{1}[\(]{4})'
        term = '.'
        self._clearStems(x, x2, term)
        y = re.compile(r'[\)]{4}([\.]{1})[\)]{1}[\.]{1}')
        y2 = r'([\)]{4})%s[\)]{1}([\.]{1})'
        term = '.'
        self._clearStems(y, y2, term)

    def getloops(self):
        self.deloopMega()
        self._loopInsert()
        startloop = [m.start() for m in re.finditer('[\.\?]{2}[\(]{4}', self.rnaStruc)]
        endloop = [n.end() for n in re.finditer('[\)]{4}[\.\?]{2}', self.rnaStruc)]
        startloopC = map (lambda x: x +3, startloop)
        endloopC = map (lambda x: x -2, endloop)
        prev = 0
        hairLoc = []
        for x,y in zip(startloopC, endloopC):
            while x > y:
                y = next(endloopC)
            if prev > x:
                hairLoc[-1][1] = y
                prev = y
                continue
            hairLoc.append([x,y])
            prev = y
        looplocs = self.gapcorrect(hairLoc)
        return looplocs, self.score

class Alignstruc(Rnatricks):
    def __init__(self, pathDir, rnafile, chunk=1, filechunk=100):
        self.align = AlignIO.read(Path(pathDir,rnafile), 'fasta')
        self.rnafile = rnafile
        self.pathDir = pathDir
        self.chunk = chunk
        self.rnaStrucAlign = MultipleSeqAlignment([])
        self.filechunk = filechunk

    def log_result(self,result):
        result_list.append(result)

    def pooladmin(self, pool, poolout = ''):
        pool.close()
        pool.join()
        pool.terminate()
        if poolout:
            pool.get()

    def poolProcess(self, poolout):
        # https://www.biostars.org/p/966/
        for kv in poolout.get():
            # rnaDict.update({''.join(kv.keys()):''.join(kv.values())})
            sr = SeqRecord(Seq(''.join(kv.values())), ''.join(kv.keys()), '', '')
            self.rnaStrucAlign.append(sr)

# https://stackoverflow.com/questions/1915170/split-a-generator-iterable-every-n-items-in-python-splitevery
    def split_every(self, n):
        i = iter(self.align)
        piece = list(islice(i, n))
        while piece:
            yield piece
            piece = list(islice(i, n))

    def rnaParallel(self, rnafold):
        # align = AlignIO.read(self.rnafile, 'fasta')
        processes_count = (cpu_count()-1)
        # alignCount = [[record, i, self.pathDir] for i, record in enumerate(self.align)]
        pool = Pool(processes_count)
        count = 1
        chunked = dict()
        chunked[count] = []
        chunkedAlign = list(self.split_every(self.filechunk))
        rnafout = Path(self.pathDir, re.sub('fa$','rnaStr.fa',str(self.rnafile)))
        if rnafout.is_file():
            rnafout.unlink()
        # pool.apply_async(gapTrans, args=(record, i, self.pathDir), callback = self.log_result)
        # poolout = pool.map_async(gapTrans, [[rC, self.rnafile] for rC in chunkedAlign], chunksize=self.chunk)
        poolout = pool.map_async(gapTrans, [[rC, rnafout, rnafold] for rC in chunkedAlign], chunksize=self.chunk)
        print ('Processing RNA secondary structures')
        rnastrucs = self.pooladmin(pool)
        rnastrucs = poolout.get()        
        hd5out = Path(re.sub('fa$','freeEnergy.hd5',str(self.rnafile)))
        if hd5out.is_file():
            hd5out.unlink()
        dd.io.save(hd5out, rnastrucs[0], compression='blosc')

# def gapTrans(record, count, path):
def gapTrans(rC):
    totalscore = []
    records = rC[0]
    rnafout = rC[1]
    rnafold = rC[2]
    for record in records:
        # /usr/local/bin/RNAfold
        seq = re.sub(r'-', '', str(record.seq))
        TUtran = {84:85}
        seq = seq.translate(TUtran)
        process = Popen(['echo %s | %s' % (seq, rnafold)],
                        stdout=PIPE, 
                        stderr=PIPE,
                        shell=True
                        )
        stdout, stderr = process.communicate()
        rnaStrucTmp = stdout.decode().split('\n')
        # rnaStruc = rnaStrucTmp[1]
        rnaStruc = rnaStrucTmp[1].split(' ')
        if len(rnaStruc) > 1: 
            if re.search(r'\(\-[0-9]+',rnaStruc[1]):
                score = float(rnaStruc[1].lstrip('(').rstrip(')'))
                totalscore.append(score)
            else:
                totalscore.append(
                    record.id + 
                    ' had problems with free energy score, please check')
        else:
            totalscore.append(record.id + ' free energy absent')
        rnatran = {46:94}
        mystruc = list(rnaStruc[0].translate(rnatran))
        location = [m.start() for m in re.finditer('-', str(record.seq))]
        for loc in location:
            mystruc.insert(loc, '-')
        with open (rnafout, 'a') as fout:
            fout.write('>{}\n{}\n'.format(record.id, ''.join(mystruc)))
    return totalscore
    
class Plotting ():
    def __init__(self, ambigList, path):
        # self.mynumpy = np.array(ambigList)
        self.ambigList = ambigList
        self.path = path

    def histogr (self, xaxis='no. ambiguous nucleotides per genome'):
        ambigNP = np.array(self.ambigList)
        # plt.hist(ambigNP, align='mid', histtype='stepfilled', bins=175)
        plt.hist(ambigNP, align='mid', histtype='stepfilled', bins=20000)
        plt.xlim([0, 200])
        # default_x_ticks = range(min(ambigNP), ambigNP[-1]+1,500)
        default_x_ticks = range(min(ambigNP), 200, 25)
        plt.title('Frequency plot (y axis) of ' + xaxis)
        # xaxis.insert(1, xaxis[0].upper())
        plt.xlabel(xaxis.capitalize())
        try:
            plt.xticks(default_x_ticks, ambigNP)
        except ValueError:
            pass  # do nothing!
        filelabel = xaxis.split()
        plt.savefig(Path(self.path, filelabel[0] + '_' + filelabel[1] +'.pdf'))
        return plt.show()

'''
if __name__ == '__main__':
    path = '/Volumes/Data/5UTR_millcuarto/testTrue1'
    path = '/Volumes/Data/'
    genFile = 'test-gen2.fasta'
    # alignFile = 'SARScov2_5UTR_1-0.75mil.aln.fa '
    ref = 'rnaref.fa'
    rnaref = Path(path, ref)
    iniQC = rnaQC(fullgenome, path)
    ambigSort, delIDs, ambigList, totalseq = iniQC.countNs()
# =============================================================================
#     
#     filehd = Path('/Volumes/Data/test-gen2.qc.hdf5')
#     myhdf5 = h5py.File(filehd, 'r')
#     item = myhdf5['5UTR']
#     # print(list(myhdf5.keys()))
#     keys = item.keys()
#     out_dict = {}
#     for key in keys:
#         out_dict[key] = item[key][()]
#     mystuff = [val.decode('utf-8') for val in out_dict.values()]
#     print (mystuff)
# =============================================================================
    iniQC = rnaQC(fullgenome, path)
    # ambigSort, delIDs, ambigList, totalseq = iniQC.countNs()
    # with open(Path(path,'genomeAmbig.pickle'), 'wb') as handle:
    #     pickle.dump(ambigList, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open(Path(path,'genomeAmbig.pickle'), 'wb') as handle:
    #     pickle.dump(ambigList, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(Path(path,'genomeAmbig.pickle'), 'rb') as handle:
        ambigList = pickle.load(handle)
    threshold = 100
    delIDs = [ambig for ambig in ambigList if ambig > threshold]
    print(len(delIDs))
    # print (len(delIDs))
    # iniPlot = Plotting(ambigList,path)
    # iniPlot.histogr()
    # iniPlot = Plotting(ambigList, path)
    # iniPlot.histogr()

    # refLoc = Path(path, ref)
    # loopIni = Rnaloop(refLoc, path)
    # rnaStruc, score = loopIni.getloops()
    # print(rnaStruc)
    
    # chunk = 2
    # rnaIni = Alignstruc(masterAlign, path, chunk)
    # rnaStruc = rnaIni.rnaParallel()
    # raise SystemExit('')
'''    
