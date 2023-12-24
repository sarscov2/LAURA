#!/usr/bin/env python3
from __future__ import print_function

# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 16:56:40 2022

@author: michaelgaunt
"""

'''
* Pull in the files
* pairwise comparison for identity ... change the structural output header to match
* not identical perform set operation to look at the differences 
* excess nucleotide versus structure

Use two output headers, the sequence header and the structure header
'''
from pathlib import Path
import re, sys, gc, csv, os, errno
from Bio import AlignIO
# from Bio import SeqIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
# import matplotlib.patches as mpatches
# import matplotlib.patches as circle
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.AlignIO.PhylipIO import SequentialPhylipWriter
from collections import defaultdict

class Settricks():
    def __init__(self, pathway, nonfitDict = '', seqFile = '', compensatoryMu = {}):
        self.pathway = pathway
        self.nonfitDict = nonfitDict
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu
    
    def trans(self):
        return {94:46}, {45:None}, {84:85}

    def stockFileN(self, cMu, count, phy=False, nonfit=False, seqgrp='', rnagrp=''):
        if phy:
            outstuff = cMu + '.S' + str(count) + '.phy'
        elif nonfit:
            outstuff = seqgrp + '_' + rnagrp + '.nonfit.' + str(count) + '.sth'            
        else:
            outstuff = cMu + '.S' + str(count) + '.sth'
        return re.sub('fa$', outstuff, Path(self.seqFile).name)
    
    def mkdirs(self, pathway):
        if not pathway.is_dir():
            os.mkdir(pathway)
            
    def _nucgaps(self,gaps,cMu,ncbi):
        # ^ to "."  T to U
        HatDottranCombi = {94:46, 84:85} 
        nogaps = {45:None}
        if str(gaps) == '1':
            return str(
                self.compensatoryMu[cMu]['nuc'][ncbi].translate(HatDottranCombi))
        else:
            return str(
                self.compensatoryMu[cMu]['nuc'][ncbi].
                translate(HatDottranCombi).translate(nogaps))

    def _mseqs(self,ncbi,gaps, cMu):
        name = ncbi
        myseq = self._nucgaps(gaps, cMu, ncbi)
        sr = SeqRecord(Seq(myseq), id=name, name=name, description='')
        return sr

    def _rnastruc(self, gaps, cMu):
        HatDottran, nogaps, _ = self.trans()
        if str(gaps) == '1':
            return str(self.compensatoryMu[cMu]['struc'].translate(HatDottran))
        else:
            return str(self.compensatoryMu[cMu]['struc'].translate(HatDottran).translate(nogaps))
        
    def _setnonStockFile(self):
        # NOTE: This file is written to via append, therefore critical 
        # identical files are removed
        pathway3 = Path(self.pathway, 'nonfitStockh')
        self.mkdirs(pathway3)
        seqFile2 = re.sub(r'qc.*.fa$', 'rnaStruc.txt', Path(self.seqFile).name)
        nonStockF = Path(pathway3, seqFile2)
        if nonStockF.is_file():
            os.unlink(nonStockF)
        return nonStockF

class Comwrit():
    def __init__(self, pathway, nonfitDict = '', seqFile = '', compensatoryMu = {}):
        self.pathway = pathway
        self.nonfitDict = nonfitDict
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu

    def _dictprep(self):
        tree2 = lambda: defaultdict(tree2)
        countseqDict = tree2()
        countrnaDict = defaultdict()
        for k in self.nonfitDict.keys():
            for k2 in self.nonfitDict[k].keys():
                for k3 in self.nonfitDict[k][k2].keys():
                    noseqs = len(self.nonfitDict[k][k2][k3]['id'])
                    if k2 in countseqDict[k]:
                        if len(countseqDict[k][k2]) > 0:
                            countseqDict[k][k2].append(noseqs)
                        else:
                            countseqDict[k][k2] = [noseqs]
                    else:
                        countseqDict[k][k2] = [noseqs]
        for k in countseqDict.keys():
            countrnaDict[k] = {}
            for k2 in countseqDict[k].keys():
                freq = sum(countseqDict[k][k2])
                countseqDict[k][k2] = {'stat': freq}
                if 'stat' in countrnaDict[k]:
                    countrnaDict[k]['stat'].append(freq)
                else:
                    countrnaDict[k]['stat'] = [freq]
        for k in countrnaDict.keys():
            freq = sum(countrnaDict[k]['stat'])
            del countrnaDict[k]['stat']
            countrnaDict[k] = {'stat':freq}
    
        countrnaDict = sorted(countrnaDict.items(), key=lambda x:x[1]['stat'], reverse=True)        

        countseqDict = [{y[0]:{'id':sorted (countseqDict[y[0]].items(), key=lambda x:x[1]['stat'], reverse=True)}} for y in countrnaDict]
        countseqDict2 = {}
        for y in countseqDict:
            for k,v in y.items():
                countseqDict2.update({k:v})
        del countseqDict
        for k in countseqDict2.keys():
            for v in countrnaDict:
                if k is v[0]:
                    countseqDict2[k].update({'tot':v[1]['stat']})
        return countseqDict2

    def writenonFit(self, masterf):
        gaps = 0
        countseqDict2 = {}
        misfit = True
        suffix = 'misFit-report.tsv'
        countseqDict2 = self._dictprep()
        fout = self._findoutfile(suffix, gaps)
        HatDottran, nogaps, TUtran = self.trans()
        tot, tot2,j,i = 0,0,0,0
        for i,k in enumerate(countseqDict2.keys()):
            i += 1
            tot += countseqDict2[k]['tot']
            if countseqDict2[k]['tot'] > 1:
                j += 1
                tot2 += countseqDict2[k]['tot']
        with open (fout, 'w') as fo:
            fieldnames, fieldnames2, csvwriter, footnote = self._writeini(
                fo, misfit, masterf)
            csvwriter.writerow(fieldnames)
            csvwriter.writerow(["", 
                                i, j, tot, tot2, tot-tot2
                                ])
            csvwriter.writerow([])
            csvwriter.writerow(fieldnames2)
            dataout = list()
            for rna in countseqDict2.keys():
                t = sum([s[1]['stat'] for s in countseqDict2[rna]['id']])
                for i,seqs in enumerate(countseqDict2[rna]['id']):
                    innerDict = self.nonfitDict[rna][seqs[0]].values()
                    ids = [a2 for a in innerDict for a2 in sorted(a['id'])]
                    sequ = ''.join(
                        [c2 for c1 in innerDict for c2 in c1['seq']]
                        ).translate(nogaps).translate(TUtran)
                    rnau = ''.join([c2 for c1 in innerDict for c2 in c1['rna']])
                    if i == 0:
                        dataout = [rna, len(countseqDict2[rna]['id']), t, seqs[0], seqs[1]['stat'], ids[0], ','.join(ids), rnau, sequ]
                    else:
                        dataout = [""] * 3 + [seqs[0], seqs[1]['stat'], ids[0], ','.join(ids), rnau, sequ]
                    csvwriter.writerow(dataout)
            for r in range(3):
                csvwriter.writerow([])
        with open (fout, 'a') as fo:                    
            fo.write(footnote)
                    
    def _findoutfile(self, suffix, gaps):
        if str(gaps) == '1':
            seqFile2 = re.sub('fa$', suffix, self.seqFile)
        else: 
            if re.search(r"\.qc\.", str(self.seqFile)):
                seqFile2 = re.sub(r'qc.*fa$', suffix, self.seqFile)
            elif re.search(r"\.dedup\.", str(self.seqFile)):
                seqFile2 = re.sub(r'dedup.*fa$',suffix, self.seqFile)
            else:
                seqFile2 = re.sub(r"fa$", suffix, self.seqFile)
        fout = Path(self.pathway, seqFile2)
        if fout.is_file():
            fout.unlink()
        return fout
    
    def _footnote(self, masterf):
        footnote = f"Please refer to '{masterf}' for master RNA structures and sequences.\n" + \
            "This file describes the 'misfit data'.\n" + \
                "These are unique RNA secondary structures compared with the master RNA secondary\n" + \
                "structure for an identical RNA sequence WITHIN THE LOCUS SPECIFIED AT INPUT \n" +\
                    "(defaults to 1-400bp otherwise) and is reported in the '...allele.fa'" +\
                        "file name.\nThis is a result of SNPs outside that defined locus (SODL) causing\n" +\
                            "changes in the RNA secondary structure within the defined locus.\n" +\
                                "Changing the total locus size of the RNA alignment to the defined" + \
                                    " locus size\nor vice versa will remove the 'misfit' entries," + \
                                        "a 'misfit' calculation is sought."
        return footnote
    
    def _writeini(self, fo, misfit=False, masterf=''):
        csvwriter = csv.writer(fo, delimiter="\t", quotechar="'")
        fieldnames = ['Summary stats', 'Total no. unique RNA strucs.', 
                      'No. unique RNA strucs >= 2 haplotypes', 
                      'Total no. seqs',
                      'No seqs >= 2 haplotypes per unique RNA struc.',
                      'No seqs with 1 haplotype per unique RNA struc.']
        if misfit:
            fieldnames2 = ['Master RNA struc. NCBI code', 
                           'Unique RNA struc. (against master struc.)',  
                           'Total seqs w. unique RNA strucs (against master)', 
                           'Master seq NCBI code', 
                           'No. NCBI seqs w. unique RNA struc. (against master)', 
                           'NCBI seq code rep w. SNPs outside defined locus (wSODL)', 
                           'All NCBI seq sSODL codes',
                           'Unique RNA Structure', 'wSODL sirus seq']
            footnote = self._footnote(masterf)
            return fieldnames, fieldnames2, csvwriter, footnote
        else:
            fieldnames2 = ['Struc. NCBI code', 'No. haplotypes', 'Total NCBI seqs', 'RNA Structure', 'NCBI seq code','No. NCBI seqs', 'Virus seq']
            return fieldnames, fieldnames2, csvwriter
    
    def writeCompMu(self):
        gaps = 0
        suffix = 'report.tsv'
        fout = self._findoutfile(suffix, gaps)
        HatDottran, nogaps, TUtran = self.trans()
        with open (fout, 'w') as fo:
            fieldnames, fieldnames2,csvwriter = self._writeini(fo)
            csvwriter.writerow(fieldnames)
            csvwriter.writerow(["", self.compensatoryMu['stats']['TotalUniqRNAstr'],
                                self.compensatoryMu['stats']['UniqRNA>=2'],
                                self.compensatoryMu['stats']['TotalSeqs'],
                                self.compensatoryMu['stats']['Seqs>=2'],
                                self.compensatoryMu['stats']['Seqs=1']
                ])
            csvwriter.writerow([])
            csvwriter.writerow(fieldnames2)
            for cMu in self.compensatoryMu:
                if cMu != 'stats':
                    if str(gaps) == '1':
                        structure = str(
                            self.compensatoryMu[cMu]['struc'].translate(HatDottran)
                            )
                    else:
                        structure = str(
                            self.compensatoryMu[cMu]['struc'].translate(HatDottran).translate(nogaps)
                            )
                    dataOut = [cMu, self.compensatoryMu[cMu]['haplotype_freq'],
                            self.compensatoryMu[cMu]['ncbi_total'],
                            structure]
                    for ncbi in self.compensatoryMu[cMu]['ncbi_freq']:
                        if str(gaps) == '1':
                            seq = self.compensatoryMu[cMu]['nuc'][ncbi].translate(TUtran)
                        else:
                            seq = self.compensatoryMu[cMu]['nuc'][ncbi].translate(nogaps).translate(TUtran)
                        dataOutAppd = dataOut + [ncbi, 
                                             self.compensatoryMu[cMu]['ncbi_freq'][ncbi], 
                                             str(seq)]
                        csvwriter.writerow(dataOutAppd)
                        dataOut = [""] * 4
        return Path(fout).name

class Comalignwrit(Settricks):
    def __init__(self, pathway, seqFile = '', compensatoryMu = {}):
        self.pathway = pathway
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu
        super().__init__(pathway, seqFile, compensatoryMu)
    
    def _phywrit(self, outAlign, phyout, cMu):
        if len(outAlign) < 1:
            return
        for record in outAlign:
            record.id = self.alignCount(record.id, record.seq, cMu)
        with open(phyout, 'w') as output_handle:  # phylipID_width
            SequentialPhylipWriter(output_handle).write_alignment(
                outAlign, id_width=35)
            
    def _MSA(self, outAlignTmp, gaps, cMu, rna):
        return MultipleSeqAlignment(outAlignTmp, column_annotations={
                    "secondary_structure": str(rna)})
    
    def _nonFitseq(self, nameList, key):
        name = "|".join(nameList)
        _, nogaps, TUtran = self.trans()
        return SeqRecord(
            Seq(key.translate(nogaps).translate(TUtran)), 
            id=name, name=name, description='')

    def alignCount(self, recordid, seqid, cMu):
        warn = ''
        genecount = str(self.compensatoryMu[cMu]['ncbi_freq'][recordid])
        if re.search(r'[NRYBDHVWSMK]', str(seqid), re.IGNORECASE):    
            warn = 'contains_N'
        return re.sub(r'_$', '', '_'.join([recordid, 'alltaxa', genecount, warn]))

class Comstock(Comalignwrit, Settricks):
    def __init__(self, pathway, nonfitDict = '', seqFile = '', compensatoryMu = {}):
        self.pathway = pathway
        self.nonfitDict = nonfitDict
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu
        super().__init__(pathway, nonfitDict, seqFile, compensatoryMu)

    def writeNonStockholm(self, nonStockF, gaps, cMu, count2, startfile, align, seqTotList, rnaTotList, masterf):
        footnote = self._footnote(masterf)
        outAlignDict = {}
        rna = self._rnastruc(gaps,cMu) 
        with open(nonStockF,'a') as f:
            if count2 == 1:
                f.write(footnote)
        for ncbi2 in self.compensatoryMu[cMu]['fullid']:
            for record in align:
                if ncbi2 == record.id:
                    tmp = (sorted(v) for v in seqTotList if (ncbi2 in v))
                    seqgrp = next(tmp)[0]
                    del tmp
                    tmp = (sorted(v) for v in rnaTotList if (seqgrp in v)) 
                    rnagrp = next(tmp)[0]
                    del tmp
                    if len(self.nonfitDict[rnagrp][seqgrp][count2]['id']) > 0:
                        self.nonfitDict[rnagrp][seqgrp][count2]['id'].append(ncbi2)
                    else:
                        self.nonfitDict[rnagrp][seqgrp][count2]['id'] = [ncbi2]
                    self.nonfitDict[rnagrp][seqgrp][count2]['rna'] = rna
                    if record.seq in outAlignDict:
                        outAlignDict[str(record.seq)].append(record.id)
                    else:
                        outAlignDict = {str(record.seq):[record.id]}

        for k in outAlignDict.keys():
            if len(self.nonfitDict[rnagrp][seqgrp][count2]['seq']) > 0:
                self.nonfitDict[rnagrp][seqgrp][count2]['seq'].append(k)
            else:
                self.nonfitDict[rnagrp][seqgrp][count2]['seq'] = [k]

        outAlignTmp2 = [self._nonFitseq(outAlignDict[k], k) for k in outAlignDict.keys()]
        seqFile3 = self.stockFileN(cMu, count2, False, True, seqgrp, rnagrp)
        outalign2 = self._MSA(outAlignTmp2, gaps, cMu, rna)
        AlignIO.write(outalign2, Path(Path(nonStockF).parent,seqFile3), 'stockholm')

    def _phyout(self):
        pathway = Path(self.seqFile).parent / 'alignments'
        self.mkdirs(pathway)
        seqFilealign = re.sub('fa$',  'align.phy', str(Path(self.seqFile).name))
        phyout = Path(pathway,seqFilealign)        
        if phyout.is_file():
            os.unlink(phyout)
        return phyout, pathway
 
    def stockholm (self, startfile, seqTotList, rnaTotList, masterf):
        gaps = 0
        count, count2 = 0, 0
        # alignment.column_annotations['secondary_structure']
        i = 0
        nonStockF = self._setnonStockFile()
        phyout, phypath = self._phyout()
        align = AlignIO.read(Path(Path(nonStockF).parents[1],startfile), "fasta")
        for cMu in self.compensatoryMu:
            if cMu != 'stats':
                count += 1
                outAlignTmp = []
                for ncbi in self.compensatoryMu[cMu]['ncbi_freq']:
                    name = ncbi
                    myseq = self._nucgaps(gaps, cMu, ncbi)
                    sr = SeqRecord(Seq(myseq), id=name, name=name, description='')
                    outAlignTmp.append(sr)
                if len(outAlignTmp) == 0:
                    count2 += 1
                    self.writeNonStockholm(nonStockF, gaps, cMu, count2, startfile, align, seqTotList, rnaTotList, masterf)
                    continue
                outAlign = MultipleSeqAlignment(outAlignTmp, column_annotations={
                    "secondary_structure": str(self._rnastruc(gaps, cMu))})
                seqFile2 = self.stockFileN(cMu, count) 
                phyFile = self.stockFileN(cMu, count, True)
                pathway2 = Path(self.pathway, 'stockholm')
                self.mkdirs(pathway2)
                sthout = Path(pathway2,seqFile2)
                AlignIO.write(outAlign, sthout, 'stockholm')
                self._phywrit(outAlign, Path(phypath, phyFile), cMu)
                i += 1 

class Comprep():
    def __init__(self, pathway, nonfitDict = '', seqFile = '', compensatoryMu = {}):
        self.pathway = pathway
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu

    def findIdent(self, idList1, idList2):
        nonidentDict = {}
        for seq, id1 in idList1.items():
            nonidentDict[id1[0]] = {'ncbiList':id1, 'onetoOne':[], 'seq':seq}
            flag = 0
            for seq2, id2 in idList2.items():
                if id1 == id2:
                    nonidentDict[id1[0]]['onetoOne'].append(True)
                else:
                   flag = 1
            if flag == 1:
                nonidentDict[id1[0]]['onetoOne'].append(False)
        return nonidentDict

    def _mkseqDict(self, inDict):
        seqstrucDict = dict()
        for seq, ids in inDict.items():
            ids.sort()
            seqstrucDict[ids[0]] = {'seq':''}
            seqstrucDict[ids[0]] = {'ncbiList':[]}
            seqstrucDict[ids[0]]['seq'] = str(seq)
            seqstrucDict[ids[0]]['ncbiList'] = ids
        return seqstrucDict

    def _combCompMuSeq (self, nucDict, strucDict):
        for k in self.compensatoryMu.keys():
            try: 
                del self.compensatoryMu[k]['bool']
                del self.compensatoryMu[k]['residual']
            except:
                pass
            # TODO: other side
            self.compensatoryMu[k].update({'haplotype_freq': len(self.compensatoryMu[k]['original'])})
            self.compensatoryMu[k].update({'nuc':{}})
            self.compensatoryMu[k].update({'nucname':{}}) 
            self.compensatoryMu[k].update({'ncbi_freq':{}})
            self.compensatoryMu[k].update({'ncbi_freq2':{}})
            self.compensatoryMu[k].update({'ncbi_total':int()})
            self.compensatoryMu[k].update({'ncbi_total2':len(strucDict[k]['ncbiList'])})
            self.compensatoryMu[k].update({'struc':strucDict[k]['seq']})
            coMuList = [st[0] for st in self.compensatoryMu[k]['original']]
            coMuLen = [len(st) for st in self.compensatoryMu[k]['original']]
            self.compensatoryMu[k]['ncbi_total'] = sum(coMuLen)
            unorderDict = dict()
            for x in zip (coMuList, coMuLen):
                unorderDict.update({x[0]:x[1]})
            orderDict = sorted(unorderDict.items(), 
                   key=lambda x: x[1], reverse = True)
            self.compensatoryMu[k]['ncbi_freq'].update(orderDict)
            cMuflat = [cMu for cMuOrig in self.compensatoryMu[k]['original'] for cMu in cMuOrig]
            for nD in nucDict.keys():
                if nD in cMuflat:
                    self.compensatoryMu[k]['nuc'].update({nD: nucDict[nD]['seq']})
                    self.compensatoryMu[k]['ncbi_freq2'].update({nD: len(nucDict[nD]['ncbiList'])})
                    self.compensatoryMu[k]['nucname'].update({nD:nucDict[nD]['ncbiList']})
        return self

    def _compyMuReorder(self):
        # TODO: Check _getHaplotypeGp
        # and change False to True
        newdict = sorted(self.compensatoryMu.items(), 
               key=lambda x: (x[1]['haplotype_freq'], x[1]['ncbi_total']), reverse = True)
        del self.compensatoryMu
        self.compensatoryMu = {nd[0]:nd[1] for nd in newdict}
        del newdict
        gc.collect()
        return self

    def removeIdent(self, seqstrucList):
        allidentList = list()
        alluniqList = list()
        for key in seqstrucList:
            if any(seqstrucList[key]['onetoOne']):
                allidentList.append(seqstrucList[key]['ncbiList'])
            else:
                alluniqList.append(seqstrucList[key]['ncbiList'])
        return allidentList, alluniqList

    def uniqCoreCalc(self, ident1,uniq1,ident2,uniq2):
        # ident1 is the query
        # ident2 is the reference
        uniqdict = dict()
        try:
            assert ident1.sort() == ident2.sort()
        except AssertionError:
            raise AssertionError('The identity check sequence data versus RNA structure failed')
            sys.exit(0)
        for uq1 in uniq1:
            uq1.sort()
            uniqdict[uq1[0]] = {'residual':[]}
            uniqdict[uq1[0]].update({'bool':[]})
            uniqdict[uq1[0]].update({'original':[]})
            uniqdict[uq1[0]].update({'fullid':uq1})
            # sequence will need to go here
            setuq1 = set(uq1)
            for uq2 in uniq2:
                uq2.sort()
                residual = setuq1 - set (uq2)
                # This simply means its empty
                if set(residual) == setuq1:
                     uniqdict[uq1[0]]['bool'].append(False)
                # The query is entirely contained in the reference
                elif residual == set():
                     uniqdict[uq1[0]]['bool'].append(None)
                else:
                # This is the surplus, compenstory mutation if
                # struc is the query
                    uniqdict[uq1[0]]['bool'].append(True)
                    uniqdict[uq1[0]]['original'].append(uq2)
                    uniqdict[uq1[0]]['residual'].append(list(residual))
        return uniqdict

class Comgraph ():
    def __init__(self, pathway, nonfitDict = '', seqFile = '', compensatoryMu = {}):
        self.pathway = pathway
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu
        # super().__init__(pathway)
# TODO: return to this section for sequence output
    def _getHaplotypeGp(self):
        # haploPerStruc = [[k, compensatoryMu[k]['original']] for k in compensatoryMu.keys()]
        haploPerStruc = [[k, self.compensatoryMu[k]['original']] for k in self.compensatoryMu.keys()]
        hapStrFreqDi = {}
        for hps in haploPerStruc:
            hapStrFreqDi[hps[0]] = {}
            noSeqs = [len(h) for h in hps[1]]
            hapStrFreqDi[hps[0]] = {'HaplotypesPerStr':len(noSeqs), 'SeqPerStr':sum(noSeqs)}
        # self._countComMu(allseqIdentLi, compensatoryMu)
        hapStrFreqDi2 = sorted(
            hapStrFreqDi.items(),
            key=lambda x: (x[1]['HaplotypesPerStr'], x[1]['SeqPerStr']), reverse = False
            )
        hapStrFreqDi3 = {hsfd2[0]:hsfd2[1] for hsfd2 in hapStrFreqDi2}
        del hapStrFreqDi, hapStrFreqDi2
        gc.collect()
        return hapStrFreqDi3

    def _mkDfDict(self, hapStrFreqDi3):
        hapStrFreqDi = {}
        StrVSeqGraph = Counter()
        for i,k in enumerate(hapStrFreqDi3.keys()):
            hapStrFreqDi[k] = {}
            if hapStrFreqDi3[k]['HaplotypesPerStr'] in StrVSeqGraph:
                StrVSeqGraph[hapStrFreqDi3[k]['HaplotypesPerStr']] += hapStrFreqDi3[k]['SeqPerStr']
            else:
                StrVSeqGraph[hapStrFreqDi3[k]['HaplotypesPerStr']] = {}
                StrVSeqGraph[hapStrFreqDi3[k]['HaplotypesPerStr']] = hapStrFreqDi3[k]['SeqPerStr']
                
            hapStrFreqDi[k] = {'HaplotypesPerStr':hapStrFreqDi3[k]['HaplotypesPerStr'],
                               'SeqPerStr':hapStrFreqDi3[k]['SeqPerStr'],
                               'StrucNo':i
                               }
        del hapStrFreqDi3
        gc.collect()
        StrVSeqGraph = sorted(
            StrVSeqGraph.items(),
            key=lambda x: x[0], reverse = False
            )
        return StrVSeqGraph, hapStrFreqDi

    def _mkXY(self, StrVSeqGraph, hapStrFreqDi):
        x2 = [svsg[0] for svsg in StrVSeqGraph]
        y2 = [svsg[1] for svsg in StrVSeqGraph]
        x = [hapStrFreqDi[k]['HaplotypesPerStr'] for k in hapStrFreqDi.keys()]
        y = [hapStrFreqDi[k]['SeqPerStr'] for k in hapStrFreqDi.keys()]
        return x2, y2, x, y

    def _dfs(self, x, y, x2, y2):
        df = pd.DataFrame({'x':x, 'y':y})
        df2 = pd.DataFrame({'x2':x2, 'y2':y2})
        # https://stackoverflow.com/questions/17709270/create-column-of-value-counts-in-pandas-dataframe
        zdf= df['x'].value_counts(
            ).sort_index().to_frame().reset_index().rename(
            columns={'index':'allelesPerStr'}
            )
        newidx = [r for r in range(df2.shape[0])]
        zdf = zdf.reindex(index=zdf.index[::-1])
        zdf['newidx'] = newidx
        zdf['SeqsPerStr'] = df2.reindex(index=df2.index[::-1])['y2']
        total = zdf['SeqsPerStr'].sum()
        xnew, ynew, ynew2 = [],[],[]
        for i, n in enumerate(zdf['allelesPerStr'].tolist()):
            zdf1 = zdf[zdf.allelesPerStr > (n-1)]
            xnew.append(n)
            ynew.append(zdf1['SeqsPerStr'].sum())
            ynew2.append(zdf1['SeqsPerStr'].sum()/total)
        dfnew = pd.DataFrame({'xnew':xnew, 'ynew':ynew, 'ynew2':ynew2})
        return df, df2, zdf, dfnew

    def _overViewG(self, df, df2, zdf):
        plt.scatter(df.x, df.y, color = '#FF0000',s=5)
        plt.yscale('log')
        plt.title('Unsummed Genbank sequences')
        plt.xlabel('No. haplotypes per RNA structure')
        plt.ylabel('Individual Genbank sequences per unique RNA structure')
        imageout = Path(self.pathway, str(self.seqFile).replace('dedup.fa', 'unsummed.pdf'))
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.savefig(imageout,format='pdf',dpi=300, pad_inches=0)
        plt.show()
        plt.close()
        plt.clf()
        self._colourPlt(zdf, df2)
        self._histoPlt(df)
        self._histoPlt2(df)
        gc.collect()

    def _numpyShape(self, dfnew):
    # TODO: Precise value 10-38 may need changing
        X = dfnew.iloc[10:38, 0].values.reshape(-1, 1)
        xr = [r for r in range(1,69)]
        dfx = pd.DataFrame({'X2':xr[::-1]})
        X2 = dfx.iloc[:, 0].values.reshape(-1,1)
        Y = dfnew.iloc[10:38, 2].values.reshape(-1, 1)
        return X, X2, Y

    def _LR(self, linear_regressor, X, X2, Y):
        linear_regressor.fit(X, Y)  
        Y_predr2 = linear_regressor.predict(X)
        Y_pred = linear_regressor.predict(X2)
        return Y_pred, Y_predr2

    def _outStats(self, linear_regressor, Y, Y_predr2):
        print(r2_score(Y, Y_predr2))
        print (str(*linear_regressor.coef_).replace('[','').replace(']',''))
        print (*linear_regressor.intercept_)
        print (str(*linear_regressor.predict([[40]])).replace('[','').replace(']',''))

    def _colourPlt(self, zdf, df2):
        my_colours = {0:'black',1:'red',2:'darkturquoise',3:'blueviolet',4:'darkorange',5:'cornflowerblue', 6:'navy'}
        colRanges = {6:[0,1], 5:[2,2], 4:[3,9], 3:[10,19],2:[20,59],1:[60,210], 0:[211,50000]}
        vocabulary = list()
        zdf = zdf.iloc[::-1]
        for struc in zdf.x.tolist():
            coltmp = [key for key, (low, high) in colRanges.items() if low <= struc <= high]
            vocabulary.append(coltmp.pop())
        plt.scatter(df2.x2, df2.y2, color = 'red', s=10)
        for x in range(0,df2.shape[0]):
            plt.scatter(df2.x2[x], df2.y2[x], 
                        color = my_colours.get(vocabulary[x], 'black'), s=15, 
                        label=my_colours.get(vocabulary[x], 'black'))
        plt.xlabel("No. of viral haplotypes per unique 5' UTR RNA structure(s)")
        plt.ylabel("No. of Genbank sequences (post-QC)")
        plt.xticks([x for x in range(0, 110, 5)], rotation = 45)
        plt.yticks([y for y in range(0, 8750, 500)])
        plt.legend(["1", "2", "3-9",  "10-19", "20-59", "60-202",  "203+"], ncol=1, title='No. RNA structure(s)')
        imageout = Path(self.pathway, str(self.seqFile).replace('dedup.fa', 'noncu.pdf'))
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.savefig(imageout,format='pdf',dpi=300, pad_inches=0)
        plt.show()
        plt.show()
        plt.close()
        plt.clf()

    def _seqVrna(self,zdf):
        plt.scatter(zdf.x, zdf.allelesPerStr, color = '#FF0000',s=5)
        plt.xscale('log')
        plt.yscale('log')
        # plt.title('Unsummed Genbank sequences')
        plt.xlabel('No. RNA structure')
        plt.ylabel('Haplotypes per unique RNA structure')
        imageout = Path(self.pathway, str(self.seqFile).replace('dedup.fa', 'strVhaps-loglog.pdf'))
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.savefig(imageout,format='pdf',dpi=300, pad_inches=0)
        plt.show()
        plt.close()
        plt.clf()
        plt.scatter(zdf.x, zdf.SeqsPerStr, color = '#FF0000',s=5)
        plt.xscale('log')
        plt.yscale('log')
        # plt.title('Unsummed Genbank sequences')
        plt.xlabel('No. RNA structures')
        plt.ylabel('Summed genbank sequences')
        imageout = Path(self.pathway, str(self.seqFile).replace('dedup.fa', 'strVseqs-loglog.pdf'))
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.savefig(imageout,format='pdf',dpi=300, pad_inches=0)
        plt.show()
        plt.close()
        plt.clf()
        gc.collect()

    def _histoPlt(self, df):
        for n in [0,1,2,3]:
            dfii = df[df.x > n]
            maxv = dfii.loc[:,'x'].max()
            plt.hist(dfii.x, bins = 105)
            if n in [1,2,3]:
                plt.xlabel("More than {} viral haplotypes per unique 5' UTR RNA structure(s)".format(n))
            else:
                plt.xlabel("All viral haplotypes per unique 5' UTR RNA structure(s)")
            plt.ylabel("No. RNA structures(s)")
            plt.xticks([x for x in range(0, 110, 5)], rotation = 45)
            imageout = Path(self.pathway, str(self.seqFile).replace('dedup.fa', 'n{}.histo.pdf'.format(n)))
            plt.yscale('log')
            plt.gcf().subplots_adjust(bottom=0.15)
            plt.savefig(imageout,format='pdf',dpi=300, pad_inches=0)
            plt.show()
            plt.close()
            plt.clf()
        del dfii, imageout
        gc.collect()

    def _histoPlt2(self, df):
        # TODO: consider the bin size
        plt.hist(df.y, bins = 212)
        plt.xlabel("Histogram of Genbank sequences (post-QC)")
        plt.ylabel("Number of RNA structures per bin (40 sequences per bin)")
        plt.xticks([x for x in range(0, 8500, 500)], rotation = 45)
        imageout = Path(self.pathway, str(self.seqFile).replace('dedup.fa', 'histo2.pdf'))
        plt.yscale('log')
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.savefig(imageout,format='pdf', dpi=300, pad_inches=0)
        plt.show()
        plt.close()
        plt.clf()
        del df
        gc.collect()

    def _plotLR(self, zdf, dfnew, X, Y, X2, Y_pred, Y_predr2):
        colRanges = {6:[0,1], 5:[2,2], 4:[3,9], 3:[10,19],2:[20,59],1:[60,210], 0:[211,50000]}
        vocabulary = list()
        for struc in zdf.x.tolist(): 
            coltmp = [key for key, (low, high) in colRanges.items() if low <= struc <= high]
            vocabulary.append(coltmp.pop())
        my_colours = {0:'black',1:'red',2:'darkturquoise',3:'blueviolet',4:'darkorange',5:'cornflowerblue', 6:'navy'}
        for x in range(0,dfnew.shape[0]):
            plt.scatter(dfnew.xnew[x], dfnew.ynew2[x], color = my_colours.get(vocabulary[x], 'black'), 
                        s=15, label=my_colours.get(vocabulary[x], 'black'))
        plt.plot(X2, Y_pred, color='grey')
        plt.xlabel("No. of viral haplotypes per unique 5' UTR RNA structure(s)")
        plt.ylabel("Population frequency (post-QC)")
        plt.xticks([x for x in range(0, 110, 5)], rotation = 45)
        plt.yticks([y/10 for y in range(1, 11)])
        plt.legend(["1", "2", "3-9",  "10-19", "20-59", "60-202",  "203+"], ncol=1, title='No. RNA structure(s)')
        imageout = Path(self.pathway, str(self.seqFile).replace('dedup.fa', 'pdf'))
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.savefig(imageout,format='pdf',dpi=300, pad_inches=0)
        plt.show()

    def plotComp(self):
        hapStrFreqDi3 = self._getHaplotypeGp()
        # hapStrFreqDi3 = self._getHaplotypeGp(compensatoryMu)
        StrVSeqGraph, hapStrFreqDi = self._mkDfDict(hapStrFreqDi3)
        x2,y2,x,y = self._mkXY(StrVSeqGraph, hapStrFreqDi)
        df, df2, zdf, dfnew = self._dfs(x, y, x2, y2)
        self._seqVrna(zdf)
        self._overViewG(df, df2, zdf)
        X, X2, Y = self._numpyShape(dfnew)
        linear_regressor = LinearRegression()
        Y_pred, Y_predr2 = self._LR(linear_regressor, X, X2, Y)
        self._outStats(linear_regressor, Y, Y_predr2)
        self._plotLR(zdf, dfnew, X, Y, X2, Y_pred, Y_predr2)

class Comuniq():
    def __init__(self, pathway, nonfitDict = '', seqFile = '', compensatoryMu = {}):
        self.pathway = pathway
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu
        # super().__init__(pathway, seqFile, compensatoryMu)        

    def getUniq (self, allseqUniqLi):
        for asUL in allseqUniqLi:
            asUL.sort()
            self.compensatoryMu[asUL[0]] = {'original':[]}    
            self.compensatoryMu[asUL[0]]['original'].append(asUL)
        return self

class Comset(Comprep, Comuniq, Comgraph, Comstock, Comwrit):
    def __init__(self, pathway, nonfitDict = '', seqFile = '', compensatoryMu={}):
        self.pathway = pathway
        self.seqFile = seqFile
        self.compensatoryMu = compensatoryMu
        tree = lambda: defaultdict(tree)
        self.nonfitDict = tree()
        super().__init__(pathway, nonfitDict, seqFile, compensatoryMu)

    def makeIDlists (self, records):
        tmprecID = dict()
        for record in records:
            tmprecID[record.seq] = []
            ids = (str(record.id)).split('|')
            ids.sort()
            tmprecID[record.seq] += ids
            if len(ids) > 1:
                record.id = '|'.join(ids)
        return tmprecID

    def _lastfile(self, fileList, searcht):
        repList = list()
        for file in fileList:
            replicate = int(str(file).split('.')[-3])
            repList.append(replicate)
        seqFinFile = [file for file in fileList if re.search(searcht + str(max(repList)) + '.dedup.fa', str(file))]
        return ''.join(str(seqFinFile[0])), replicate

    # FIXME: this is detached
    def _countComMu(self, allseqIdentLi, identNonidentLi):
        fulllist = [identNonidentLi[k]['original'] for k in identNonidentLi.keys()]
        # print (len([compensatoryMu[k]['original'] for k in compensatoryMu.keys()]))
        fulllist2 = [f for fl in fulllist for f in fl]
        fulllist3 = [i2 for f2 in fulllist2 for i2 in f2]
        print (len(fulllist), "whats this")
        print(len(fulllist2), "and this")
        print (len(fulllist3), "here too")
        ident1 = [asi for asi in allseqIdentLi]
        ident2 = [idi for asiLi in allseqIdentLi for idi in asiLi]
        print (len(ident1))
        print (len(ident2))

    def openoutput(self):
        fL = Path(self.pathway).glob(r'*[0-9][0-9]*.dedup.fa')
        rL = Path(self.pathway).glob(r'*[0-9][0-9]*.dedup.fa')
        fileList = [file for file in fL if not re.search('rnaSt', str(file))]
        try:
            seqFile, num = self._lastfile(fileList, '')
        except Exception as error:
            print('There are no fasta file named \'{}.dedup.fa\'.\nTechnical report: {}'.format(num,repr(error)))
            sys.exit(0)
        print ("Found LAURA sequence file {}".format(Path(seqFile).name))        
        rL2 = Path(self.pathway).glob(r'*rnaStr.[0-9]+.dedup.fa')
        rnaList = [file for file in rL2]
        if len(rnaList) == 0:
            rnaList = [file for file in rL if re.search('rnaStr', str(file))]# if re.search('rnaSt', str(file))]
        try: 
            rnaFile, _ = self._lastfile(rnaList, 'rnaStr.')
        except Exception as error:
            print('There are no structure file named \'rnaStr.{}.dedup.fa\'.\nTechnical report: {}'.format('x',repr(error)))
            sys.exit(0)
        print ("Using LAURA structure file {}".format(Path(rnaFile).name))
        if not Path(seqFile).exists():
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 
                str(seqFile))
        if not Path(seqFile).exists():
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 
                str(rnaFile))
        rnaObj =  AlignIO.read(Path(self.pathway,rnaFile),'fasta')
        seqObj =  AlignIO.read(Path(self.pathway,seqFile),'fasta')
        return seqObj, rnaObj, seqFile, rnaFile

    def _basicstats(self, allseqUniqLi, allseqIdentLi):
        a = len([len(
            self.compensatoryMu[k]['original']) for k in self.compensatoryMu.keys()]
            )
        print (
            "Total no. unique RNA structures, {}".format(
                a))
        b = a - len(allseqIdentLi)
        print ("No. unique RNA structures with >= 2 haplotypes, {}".format(b))
        c = len([asUL for uniqRNA in allseqUniqLi for asUL in uniqRNA])
        d = len([asI for asIl in allseqIdentLi for asI in asIl])
        print (
            "Total no. seqs {}".format(c+d))
        print (
            "No. seqs with >= 2 haplotypes per unique RNA struct., {}".format(
                c))
        print (
            "No. seqs with 1 haplotypes per unique RNA struct., {}".format(
                d))
        return a, b, c, d

    def lauraCom (self, startfile):
        seqObj, rnaObj, seqFile, rnaFile = self.openoutput()
        self.seqFile = seqFile
        # Find indentity between 
        # sequence fasta file and structure fasta file
        seqidDict = self.makeIDlists(seqObj)
        rnaStrIdDict = self.makeIDlists(rnaObj)
        nucDict = self._mkseqDict(seqidDict)
        strucDict = self._mkseqDict(rnaStrIdDict)
        del seqObj, rnaObj
        # Split into lists of non-identity and identity
        # for the structure - sequence comparison:
        # * part 1 is to generate the feeder files
        nonidentSeqDict = self.findIdent(seqidDict,rnaStrIdDict)
        allidentStrList = self.findIdent(rnaStrIdDict, seqidDict)
        seqTotList = tuple(v for k,v in seqidDict.items())
        rnaTotList = tuple(v for k,v in rnaStrIdDict.items())
        del seqidDict, rnaStrIdDict
        # * part 2 is to obtain the sseparate lists
        allseqIdentLi, allseqUniqLi = self.removeIdent(nonidentSeqDict)
        allStruIdentLi, allrnaStrUniqLi = self.removeIdent(allidentStrList)

        # obtain compensatory mutations against 'faulty mutations'
        faultyMu = self.uniqCoreCalc(
            allseqIdentLi,allseqUniqLi,allStruIdentLi,allrnaStrUniqLi
            )
        self.compensatoryMu = self.uniqCoreCalc(
            allStruIdentLi,allrnaStrUniqLi,allseqIdentLi,allseqUniqLi
            )
        self.getUniq(allseqIdentLi)
        del nonidentSeqDict, allidentStrList, allStruIdentLi, allrnaStrUniqLi 
        gc.collect()
        # self.plotComp(compensatoryMu)
        self._combCompMuSeq(nucDict, strucDict)

        a, b, c, d = self._basicstats(allseqUniqLi, allseqIdentLi)
        del allseqUniqLi, allseqIdentLi

        self._compyMuReorder()
        self.compensatoryMu.update({'stats':{}})

        self.compensatoryMu['stats'].update({
            'TotalUniqRNAstr':a, 'UniqRNA>=2':b, 'TotalSeqs':c+d,'Seqs>=2':c, 'Seqs=1':d
            })

        masterf = self.writeCompMu()
        self.stockholm(startfile, seqTotList, rnaTotList, masterf)
        self.writenonFit(masterf)
        return seqFile, rnaFile 

if __name__ == '__main__':
    pathway = '/Volumes/Data/5UTR_millcuarto/fullTest/fullTest2/'
    retain = True
    startfile = 'SARScov2_5UTR_1-0.75mil.full.aln2.qc.0.fa'
    iniC = Comset(pathway)
    iniC.lauraCom(startfile)
    
    