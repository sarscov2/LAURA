#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Mar 15 11:26:07 2020
TaxaSwitch is the file name single because it is the first class in the file.
It performs two different functions
1. There is the base class
2. Converts a phylip file to fasta format and vice versa
3. Slice a section of an aligment and check it for absolute identity.
   Taxon names that are identical are stacked into the taxon ID
4. Take a list of taxa of a given order and impose that ordered list onto a phylip file.
   It should complain if there is a taxon in the list but not in the alignment

@author: Michael Gaunt

"""

from __future__ import print_function
from Bio import AlignIO
# change this to from Bio.SeqIO import parse
from Bio.AlignIO.PhylipIO import SequentialPhylipWriter
import re, json, subprocess, os, gc, io, mmap
from pathlib import Path
#from abc import ABCMeta, abstractmethod
from collections import defaultdict
import Bio.Align
from datetime import datetime
from numpy import base_repr
from itertools import chain
# from shutil import copyfile
from shutil import copy2
import pandas as pd

# class AbstractClass (metaclass=ABCMeta):
#     # Class Taxatricks
#     @abstractmethod
#     def fileDel(self, outfile):
#         raise NotImplementedError
#     @abstractmethod
#     def _inOut(self, inOut, fileType):
#         raise NotImplementedError
#     @abstractmethod
#     def readFasta(self):
#         raise NotImplementedError
#     @abstractmethod
#     def readPhylip (self):
#         raise NotImplementedError
#     @abstractmethod
#     def openit (self):
#         raise NotImplementedError
#     @abstractmethod
#     def writealign(self):
#         raise NotImplementedError
#     @abstractmethod
#     def seqDedupFasta(self, fileIn, dedup_records):
#         raise NotImplementedError
#     @abstractmethod
#     def alleleWrite (self, pathOut, dedup_records):
#         raise NotImplementedError
#     # Class Convert
#     @abstractmethod
#     def phy2fa(self):
#         raise NotImplementedError
#     @abstractmethod
#     def fa2phy(self):
#         raise NotImplementedError
#     # Class Deduplication
#     def alignSlice(self, loc):
#         raise NotImplementedError
#     @abstractmethod
#     def seqDedup(self, loc):
#         raise NotImplementedError
#     # Taxaswitch
#     @abstractmethod
#     def _masterfile(self):
#         raise NotImplementedError
#     @abstractmethod
#     def seqMatch(self, loc):
#         raise NotImplementedError
#     # Class Load
#     @abstractmethod
#     def files(cls, locPhylip, locFasta, pathDir):
#         raise NotImplementedError
class TaxaBase ():
    pass

class Taxatricks ():
    def __init__(self):
        self.locPhylip = ''
        self.pangotype = ''
        self.pathDir = ''
        self.inputName = ''
        self.noiseflag = 0
        self.align = {}
        self.IDdict = {}
        self.jsonlist = {}
        self.pangoFile = ''
        self.phyalignT = {}

    def fileDel(self, outfile, *args):
        if outfile.is_file():
            outfile.unlink()
            if int(self.noiseflag):
                print('%s present, will be deleted and rewritten' %
                  (Path(outfile).name))
        else:
            return

    def _pathsort(self):
        # Just want to check whether the given directory is real otherwise stop
        if self.pathDir != '':
            tmp = self.pathDir
            self.pathDir = Path(self.pathDir)
            if self.pathDir.is_dir():
                return self.pathDir
            else:
                raise Exception('Input directory structure %s is not valid' % (tmp))
    # The following section is simply to remove an annoying bug in using Dropbox as a automated syncing system
    # There are perhaps better ways to explicitly search for Dropbox in the pathname
        else:
            mycwd = Path.home().joinpath(Path.cwd().parent.name, Path.cwd().name)
            # the method name drops the file directory whilst parent name drops the parent directory of cwd
            if mycwd.is_dir():
                return mycwd
            else:
                # If it is a hassle just give up and use the current working directory
                return Path.cwd()

# The _inOut module is for switching suffixes on out files and simply connects the file name with the directory name
# for the infile file. It works for phylip, fasta and text files using pathlib. It connects to _pathsort which is used
# to overcome an issue with Dropbox syncing and check the directory path name
    def _inOut(self, fileType, inOut, textfile, *args):
        pathDirT = self._pathsort()
        # Here we use the phylip or fasta input file for the first json output name, fasta or phylip variables will be blank hence the comprehension
        # self.locJson then goes into filedict
        if textfile == 'no':
            locJson = 'na'
        else:
            # The any statement is to make sure it is not locFasta,
            # as a pangotype or empty file for phylip
            locJson = re.sub(r'\.[a-z|A-Z]{2,4}$',r'.json',''.join([x for x in [str(self.pangotype), str(self.locPhylip)] if x != '1' and x != '0' and x != '']))
        filedict = {'phy': self.locPhylip, 'fa': self.pangotype, 'json': locJson}
        for fileformat in filedict.keys():
            if re.search(fileType, fileformat):
                # INPUT files to get the infile name and path
                if re.search(r'[I|i]n', inOut):
                    fileIn = Path.joinpath(pathDirT, filedict[fileformat])
                    if Path(fileIn).is_file():
                        if int(self.noiseflag):
                            print('Found %s' % (Path(fileIn).name))
                        return fileIn
                    else:
                        raise Exception('I cannot find the input file %s. Please check the path name.' % (fileIn))
                # OUTPUT files to get outfile name and path
                elif re.search(r'[O|o]ut', inOut):
                    if fileformat != 'json':
                        # the code below will switch the phy suffix for the fas suffix
                        # opposite = fileformat.translate(str.maketrans('phy','fas'))
                        opposite = ''.join(x for x in ['fa', 'phy'] if x != fileformat)
                        for y in args:
                            if y == 'keep':
                                opposite = fileformat
                        # connect the path of the input file with a modified output file
                        # this includes the word 'slice' and the coordinates of the slice, whilst changing the suffix as above
                        # if re.search('.dedup.', filedict[fileformat]):
                        #     filedict[fileformat] = re.sub(r'[0-9]*.dedup.','',filedict[fileformat])
                        outfileTrue = Path.joinpath(pathDirT, re.sub('.' + fileformat, self.inputName, filedict[fileformat]) + opposite)
                        self.fileDel(outfileTrue)
                        return outfileTrue
                    else:
                        try:
                            # Return and outfile if it is a text file
                            outfileTrue = Path.joinpath(pathDirT, re.sub('.' + fileformat, self.inputName, filedict[fileformat]) + 'json')
                            return outfileTrue
                        except Exception:
                            raise Exception ("Formating problems"
                                             "with the out file")
                else:
                    raise Exception(
                        "Kindly check whether you're reading in or out, %s"
                        % (inOut))
            # The elif ensures if on the first pass the filetype
            # and the dictionary variable don't match try again
            elif not re.search(fileType, fileformat):
                continue
            else:
                raise Exception('Error in Taxatricks, check valid file type'
                                ' or in/out format')

    def readFasta(self):
        # The final variable in the _inOut input is whether
        # it is a standard text file such cvs
        if self.pangotype == '':
            raise Exception('Fasta format input was omitted. \n'
                            'If %s is the file it was input'
                            'in the wrong place.'
                            % (self.locPhylip))
        fileIn = self._inOut('fa', 'in', 'no')
        try:
            return AlignIO.read(fileIn, "fasta")
        except:
            raise Exception('Fasta input file error with %s' % (Path(fileIn).name))

    def readPhylip(self):
        # The final variable in the _inOut input is whether it is a standard text file such cvs
        if self.locPhylip == '':
            raise Exception('Phylip format input was omitted. \nIf %s is the file it was input in the wrong place.' % (self.pangotype))
        fileIn = self._inOut('phy', 'in', 'no')
        try:
            return AlignIO.read(fileIn, "phylip-relaxed")
        except:
            raise Exception('Phylip input file error')
 
    def openit(self, fileIn):
        # The seqMatch method uses Fin or FileIn from
        # seqDedup, therefore _inOut is bypassed
        # fileIn = self._inOut('txt','in', 'yes')
        try:
            if fileIn.is_file():
            # TODO: consider reading from pickle
                with open(fileIn, 'rt') as t:
                    namedFile = t.readlines()
                return namedFile
        except Exception:
            raise Exception('%s was not found within class Taxaswitch' % (fileIn))

# TODO: needs to be removed - not used
    def writealign(self, outFormat, inFormat, *args):
        Fout = self._inOut(inFormat, 'out', 'no', *args)
        if re.search(r'fa', outFormat):
            AlignIO.write(self.align, Fout, outFormat)
            return Fout
        # writing out in relaxed phylip format is not trivial, requiring
        # SequentialPhylipWriter
        # see the legendary write_trees subroutine/function
        # within RYalign_AMWS.py for details
        elif re.search('phy', outFormat):
            with open(Fout, 'w') as output_handle:
                SequentialPhylipWriter(output_handle).write_alignment(self.align, id_width=50)
                return Fout
        elif re.search('nex', outFormat):
            output_handle = open(Fout, "w")
            AlignIO.write(self.align, output_handle, "nexus", "protein")
            output_handle.close()
            return Fout
        else:
            raise Exception('%s is not recognised' % (outFormat))

    def seqDedupFasta(self, fileIn, dedup_records):
        for record in fileIn:
            # Central routine using the sequence
            # as the key and append the taxa name
            dedup_records[str(record.seq)].append(record.id)
        del fileIn
        gc.collect()
        return dedup_records

    def pangoCheck (self, pangotime):
        if not pangotime.isnumeric():
            pangotime = 0
            print ('Pangotime is now 0.')
        # This is a nice bit of code pangotime has three inputs
        # if there is a value outside this it gets set to 0
        if 0 <= int(pangotime) <= 2:
            pass
        else:
            pangotime = 0
        return pangotime

    def ncbiHeadCode(self, IDdictKey, x):
        # I particularly like this way of unfolding 
        # base36 compresion however fastawrite 
        # is simpler without counting
        decoded = [self.IDdict[ID] for ID in IDdictKey]
        freqTaxa = [x.split(ID) for ID in decoded]
        flattened = list(chain.from_iterable(freqTaxa))
        return flattened
    
    def _pangoDict (self):
        z = re.compile(r'_')
        pangodicti = {}
        if Path(self.pangoFile).is_file:
    # TODO: consider reading from pickle
            with open (self.pangoFile, 'r') as din:
                for line in din:
                    if line[0] == '>':
                        line = line.lstrip('>').rstrip('\n')
                        zList = z.split(line)
                        if len(list(zList[0].split('.'))) > 1:
                            ncbiKey = re.sub(r'\.[0-9]*$', '', zList.pop(0))
                        else:
                            ncbiKey = zList.pop(0)
                        ptype = ''.join(zList[0])
                        ptypeCr = re.sub(r'([A-Z]+\.[0-9]+\.?[0-9]*)\.?[0-9]*$',r'\1', ptype)
                        pangodicti[ncbiKey] = {'pangotype':[ptype, ptypeCr]}
                        if len(zList) == 1:
                            continue
                        else:
                            country = ''.join(zList[1])
                            ptime = zList[2]
                            species = zList[3]
                            pangodicti[ncbiKey].update({'country':country})
                            pangodicti[ncbiKey].update({'ptime':ptime})
                            pangodicti[ncbiKey].update({'species':species})
            # Development
            # with open(Path(Path(self.pathDir).parent, 'truth.json'), 'a') as f:
            #     json.dump(pangodicti, f)
            return pangodicti
        else:
            raise SystemExit ('Pangotime options 1 and 2 cannot find original fasta file.')

    def pangoProc (self, pangotime, jsonlist, IDsList, ncbiKey, y, pangodicti):
        # Note flattened is the same as IDsList
        # The output is the JSON list that feeds 
        # into the histogram
        pangotypeTrun,pangotype, country, ptime, species = [],[],[],[],[]
        for IDs in IDsList:
            pangotypeTrun.append(pangodicti[IDs]['pangotype'][1])
            pangotype.append(pangodicti[IDs]['pangotype'][0])
            if str(pangotime) == '2':
                    country.append(pangodicti[IDs]['country'])
                    ptime.append(pangodicti[IDs]['ptime'])
                    species.append(pangodicti[IDs]['species'])
        jsonlist[ncbiKey] = {'pangotypeTrun':pangotypeTrun}
        jsonlist[ncbiKey].update({'pangotype':pangotype})
        jsonlist[ncbiKey].update({'freq':len(pangotype)})     
        if str(pangotime) == '2':
            jsonlist[ncbiKey].update({'country':country})
            jsonlist[ncbiKey].update({'ptime':ptime})
            jsonlist[ncbiKey].update({'species':species})
        return jsonlist

    def pdHistogram (self, epiList, istime):
        if int(istime):
            ab = re.compile(r"^([^-]+?-[^-]+?)-.*$")
            epiList = [ab.sub(r"\1", date) for date in epiList]
        epiFreqTmp = []
        df = pd.Series(epiList)
        for tup in df.value_counts().items():
            epiFreqTmp.append(str(tup[0]) + r"=" + str(tup[1]) + r";")
        epiFreq = ''.join(epiFreqTmp)
        epiFreq = re.sub(r";$", "|", epiFreq)
        return epiFreq

    def alleleWrite(self, pathOut, dedup_records, pangotime):
        pathOut = re.sub(r"dedup_tmp.fa", r"pango_v" + pangotime + "." + "allele.fa", str(pathOut))
        pathOut2 = re.sub(r"[\w]{0,4}\.[0-9]+\.dedup[^.]*\.", "", str(pathOut))
        x = re.compile('\|')
        y = re.compile(r'_')
        jsonlist = {}
        if Path(pathOut).is_file():
            self.fileDel(pathOut)
        pangotime = self.pangoCheck(pangotime)
        if str(pangotime) == '1' or str(pangotime) == '2':    
            pangodicti = self._pangoDict()
        with open(pathOut2, 'w') as output:
            for sequence, ids in dedup_records.items():
                if not self.development:
                    ids = ''.join(ids)
                IDdictKey = x.split(ids)
                flattened = self.ncbiHeadCode(IDdictKey, x)
                ncbiKey = flattened[0]
                if str(pangotime) == '0':
                    output.write(">{}_alltaxa:{}\n{}\n".format(flattened[0], len(flattened), sequence))
                if str(pangotime) == '1' or str(pangotime) == '2':
                    jsonlist = self.pangoProc(pangotime, jsonlist, flattened, ncbiKey, y, pangodicti)
                    pangofreq = self.pdHistogram(jsonlist[ncbiKey]['pangotype'], 0)
                    pangofreqTrun = self.pdHistogram(jsonlist[ncbiKey]['pangotypeTrun'], 0)
                    if str(pangotime) == '1':
                        output.write(">{}|alltaxa={}|pangotype-freq-trunc:{}pangotype-freq:{}\n{}\n".format(ncbiKey, len(flattened), pangofreqTrun, pangofreq, sequence))
                    if str(pangotime) == '2':
                        country = self.pdHistogram(jsonlist[ncbiKey]['country'],0)
                        ptime = self.pdHistogram(jsonlist[ncbiKey]['ptime'],1)
                        # species should be used when the final data 
                        # set is downloaded
                        # species = self.pdHistogram(jsonlist[ncbiKey]['species'],0)
                        output.write(">{}|alltaxa={}|pangotype-freq-trunc:{}pangotype-freq:{}country:{}time:{}\n{}\n".format(ncbiKey, len(flattened), pangofreqTrun, pangofreq, country, ptime, sequence))

    def seqwrite(self, pathOut, dedup_records):
        self.fileDel(pathOut)
    # TODO: write to pickle
        with open(pathOut, 'w') as output:
            for sequence, ids in dedup_records.items():
                output.write(">{}\n".format('|'.join(ids))+sequence + '\n')
        return dedup_records

    #  _inOut is needed here to obtain Fout, thats why it is there
    # name and seq can't be replaced however because its a loop
    def fastawrite(self, Fout, name, seq, verticalBar, x):
        name2 = ''
        identIDs = []
        strSeq = ''.join(seq)
        identIDs = x.split(name)
        if verticalBar == 1:
            name2 = '|'.join([self.IDdict[IDs] for IDs in identIDs])
        elif verticalBar == 0:
            name2 = self.IDdict[name]
        else:
            raise Exception('There is a problem with %s' % (verticalBar))
    # TODO: write to pickle
        with open(Fout, 'a') as p:
            p.write('>' + name2 + '\n' + strSeq + '\n')
            name = ''
        # returning the first sequence ID back to the calling method
                # decoded = [self.IDdict[ID] for ID in IDdictKey]
                # freqTaxa = [x.split(ID) for ID in decoded]
                # flattened = list(chain.from_iterable(freqTaxa))

    def writejson(self, jsonlist):
        Fout = self._inOut('json', 'out', 'text')
        with open(Fout, 'a') as f:
            json.dump(jsonlist, f)

    def readjson(self):
        Fin = self._inOut('json', 'in', 'text')
        with open(Fin) as f:
            idcodedict = json.load(f)
        return idcodedict

# class Convert(Taxatricks):
#     def __init__(self, locPhylip, locFasta, pathDir, IDdict, phyalignT, pangoFile, align=''):
#         # TODO: check the super method
#         # locFasta is not represented as pangotype here
#         super().__init__()
#         self.locPhylip = locPhylip
#         self.locFasta = locFasta
#         self.pathDir = pathDir
#         self.align = align
#         # The variable below is activately used in alignSlice and required by _inOut for the outfile. It is not required for
#         # for straight file conversions and therefore
#         self.inputName = ""
#         if str(locFasta) == '0' or str(locFasta) == '1':
#             self.noiseflag = locFasta
#         else:
#             self.noiseflag = 1

#     def fa2phy(self):
#         # Writing to relaxed phylip requires SequentialPhylipWrite and there is not other way
#         # It is therefore separate from alignSlice method and should remain distinct to avoid confusion
#         # Reading the fasta file is a standard Biopython method
#         self.align = super().readFasta()
#         Fout = super().writealign('phylip', 'fa')
#         if int(self.noiseflag):
#             print(f"{'Fasta to relaxed Phylip conversion is complete for %s' % (Path(Fout).name)}")
#         return Fout

#     def phy2fa(self):
#         # This is identical to alignSlice except the slicing step is omitted
#         # It has been separately written to avoid confusion and ensure the Convert class
#         # of fasta to phylip, fasta to nexus (when written), phylip to fasta etc ... is a coherant class
#         self.align = super().readPhylip()
#         super().writealign('fasta', 'phy')
#         return f"{'Phylip to fasta format conversion complete'}"

#     def getseqphylip(self, seq):
#         self.align = super().readPhylip()
#         return ''.join([str(record.seq) for record in self.align if record.id == seq])

class Deduplicate(Taxatricks):
    def __init__(self, locPhylip, noiseflag, pathDir, IDdict, phyalignT, pangoFile, alignTmp='', align='', directInput=0):
        # TODO: check the way alignTmp, align and directInput
        super().__init__(locPhylip, noiseflag, pathDir, IDdict, phyalignT, pangoFile)
        self.locPhylip = locPhylip
        self.noiseflag = noiseflag
        self.pathDir = pathDir
        self.phyalignT = phyalignT
        self.alignTmp = alignTmp
        self.align = align
        self.directInput = directInput
        self.pangoFile = pangoFile
        self.development = 0

    # if Load.files(variables).alignSlice is called it will trigger __repr__
    def alignSlice(self, loc):
        # The map function switches the string input into integers , in Python 3 map must be inside the list method
        locList = list(map(int, re.split(r'-', loc)))
        alignTmp = self.phyalignT
        if self.directInput == 1:
            fileT = ''.join([x for x in [self.locPhylip, self.locFasta] if x != ''])
    # TODO: pickle read in should go here
            if re.search(r'.phy$', fileT):
                alignTmp = super().readPhylip()
            elif re.search ('.fa', fileT):
                alignTmp = super().readFasta()
            else:
                raise Exception ('Can not determine file formate for align slicing, file must be either .phy or .fa')
        # This is the central code in the class, the rest is boiler plate
        align = alignTmp[:, locList[0]:locList[1]]
        del alignTmp
        # self.inputName = '.slice.locus' + str(locList[0]) + '-' + str(locList[1]) + "."
        # super().writealign('fasta', 'phy') # Note 'fasta' is the output and 'phylip is the input'
        if int(self.noiseflag):
            print('Slicing complete')
        return align

    # @property and its setter decorators define the sanity limits on the loc value input
    @property
    def loc(self):
        return self._loc

    @loc.setter
    def loc(self, locvalue):
        if not re.search(r'[0-9]-[0-9]', locvalue) and int(re.sub(r'-[0-9]+', '', locvalue)) < int(re.sub(r'[0-9]+-', '', locvalue)):
            raise ValueError('There is a problem with the alignment range input')
        else:
            self._loc = locvalue

    def seqDedup(self, loc, alleleflag, pangotime):
        # Using method chaining to automate the output of slicing into this method
        # Slicing is therefore a method that does not need to be calle
        align = self.alignSlice(loc)
        self.inputName = '.locus' + loc + '.dedup_tmp.'
        # Npte although the input is fasta due to the switch format option must be 'phy'
        pathOut = super()._inOut('phy', 'out', 'no')
        # initiate an approximate dict container not require importing method
        dedup_records = defaultdict(list)
        # This is the central method within Taxatricks. The readin method using SeqIO is justified
        # for Class Taxatricks, the use of a sequence as the key and taxa name as the value should be
        # here but as part of a loop it wouldn't be efficient to separate the two
        dedup_records = super().seqDedupFasta(align, dedup_records)
        # This is the write out loop for a SeqIO object, AlignIO will not work here
        if int(self.noiseflag):
            print('Temporary output for sequence deduplication complete')
        if int(alleleflag):
            super().alleleWrite(pathOut, dedup_records, pangotime)
            return dedup_records
        if self.development:
            super().seqwrite(pathOut, dedup_records)
            return pathOut
        else:
            return dedup_records

    def __repr__(self):
        return f"{self.__class__.__name__} (phyalignT={self.phyalignT!r}, phylip={self.locPhylip!r}, fasta={self.locFasta!r}, path1={self.pathDir!r} )"
    def __str__(self):
        return f"{self.phyalignT} \n{self.locPhylip} \n{self.locFasta} \n{self.pathDir} "

class Taxaswitch(Deduplicate, Taxatricks):
    def __init__(self, locPhylip, noiseflag, pathDir, IDdict, phyalignT, pangoFile):
        super().__init__(self, locPhylip, noiseflag, pathDir, IDdict, phyalignT, pangoFile)
        self.locPhylip = locPhylip
        self.pathDir = pathDir
        self.IDdict = IDdict
        self.phyalignT = phyalignT
        self.inputName = ''
        self.noiseflag = noiseflag
        self.pangoFile = pangoFile
        self.development = 0
        
    def _mkIDlist(self, name, x):
        IDname = x.split(name)
        IDsList = []
        #This section reconnects the internal label with
        # the dictionary
# TODO: Check the JSON base36 decompression with seqwrite
        ##################################
        for IDs in IDname:
            # if not int(pangotime):
            innerList = x.split(self.IDdict[IDs])
            checkNoBar = innerList[0]
            # pangotime good to here
            if len(innerList) > 1:
                for inner in innerList:
                    if x.search(inner):
                        IDname = x.split(inner)
                        for inner2 in IDname:
                            IDsList.append(inner2)
                    IDsList.append(inner)
            elif x.search(checkNoBar):
                innerList3 = x.split(checkNoBar)
                for inner3 in innerList3:
                    IDsList.append(inner3)
            else:
                IDsList.append(innerList[0])
        ##################################
        return IDsList

    def procJson (self, name, seq, loc, x, pangotime):
        # TODO: pangotime goes here
        pangotime = super().pangoCheck(pangotime)
        jsonlist = {}
        strSeq = ''.join(seq)
        # This line could be used within Taxatricks directly because
        # it is duplicated alignSlice
        locList = list(map(int, re.split(r'-', loc)))
        allele = ''
        allele = strSeq[locList[0]:locList[1]]
        IDsList = self._mkIDlist(name, x)
        y = re.compile('_')
        try:
            if not int(pangotime):
                jsonlist[IDsList[0]] = {'idlist':IDsList}
                # jsonlist[self.IDdict[identIDs[0]]].update({'seq':strSeq})
                jsonlist[IDsList[0]].update({'freq':len(IDsList)})
                # jsonlist[IDsList[0]].update({'allele':allele})
            # Note this means pangotime 1 or 2
            if int(pangotime):
                jsonlist = super().pangoProc(pangotime, jsonlist, IDsList, y)
        except:
            raise Exception ('There is a problem with %s in the procJson method' % (jsonlist))
        super().writejson(jsonlist)

    def _masterfile(self, Fin):
        master_name = []
        if self.development:
    # TODO: to check to compatibility with pickle opening
            name_file = super().openit(Fin)
            y = re.compile('>')
            for line2 in name_file:
                if y.match(line2):
                    master_name.append(line2.lstrip('>').rstrip('\n'))
                else:
                    continue
        else:
            name_file = Fin
            for fastaID in name_file.values():
                master_name.append(''.join(fastaID))
        del Fin
        gc.collect()
        return master_name

    # if set(self.idcodedict.keys()) == set(idcodecheck.keys()):
    def _seqMatchIter(self, align, master_name, Fout, loc, jsonflag, pangotime):
        super().fileDel(Fout)
        i = 0
        # Need to generate two output files, one essential contain the sequence and the other with
        # Full hearder output. These need to be deleted before output commences because they
        # are subject to an append operation
        # Fout2 = Path(re.sub(r'\.fa$', '.idTrun.fa', str(Fout)))
        Path(Fout).touch()
        x = re.compile('\|')
        # A dictionary self.idcodedict will keep count of the number of taxa sharing an identical sequence for the region specified
        # The dictionary name is just the first sequence in the duplicate fasta file
        # A dictionary is required to check all the taxa names not included in the above dictionary
        # This list is then screened against the complete taxa ids to make sure all taxa are included
        # In other words none have gone missing in the dictionary
        # idcodecheck = {}
        for record in self.phyalignT:
            for name in master_name:
            # It is critical to remove > and to a lesser extent \n tags for fasta for searching and id names. These
            # are added back later to control their position
                if self.development:
                    name = name.lstrip('>').rstrip('\n')
                name2 = ''
                if x.search(name):
                    # (?:\|), but the old style certain works
                    taxalist = x.split(name)
                    name2 = taxalist[0]
                    verticalBar = 1
#                    name2 = re.sub(y,r'\1', name, flags=re.IGNORECASE)
                    if record.id == name2:
                        i += 1
                        super().fastawrite(Fout, name, record.seq, verticalBar, x)
                        if int(jsonflag):
                            self.procJson(name, record.seq, loc, x, pangotime)
                    else:
                        # idcodecheck[name2] = 'multiple taxa check'
                        continue
                else:
                    if name == record.id:
                        # If the taxa ID lacks a '|' there is only one such sequence for the region under analysis
                        verticalBar = 0
                        i += 1
                        super().fastawrite(Fout, name, record.seq, verticalBar, x)
                        if int(jsonflag):
                            self.procJson(name, record.seq, loc, x, pangotime)
                    else:
                        # idcodecheck[name] = 'single taxa check'
                        continue
        if int(self.noiseflag):
            print ('Found %s unique sequences' % (i))
        return

    def seqMatch(self, loc, jsonflag, alleleflag, pangotime, *args):
        Fout, master_name = '', []
        # The gotofile is an option for using an exisiting dedup_tmp file to create dedup_tmp
        # This may not be needed in the final version, but was used in development
        # Parallelisation would also help.
        if str(pangotime) == '0':
            self.pangoFile = ''
        # gotofile = ''.join([z for z in args if re.search(r'dedup_tmp',z)] )
        # if gotofile:
        # # TODO: to check to compatibility with pickle opening
        #     master_name = super().openit(Path(self.pathDir, gotofile))
        #     if int(self.noiseflag):
        #         print ('Found the temporary deduplication file.')
        # else:
        #This is the standard output
        if int(self.noiseflag):
            print ('Initiating temporary deduplication file, including \'dedup_tmp\'.')
        Fin = super().seqDedup(loc, alleleflag, pangotime)
        self.inputName = '.locus' + loc + ".dedup."
        Fout = self._inOut('phy', 'out', 'no')
        # An additional complication arises to create the master file for dedup_tmp
        # if this has not been initially supplied. This needs to be tested
        # if not re.search('dedup_tmp', gotofile):
        master_name = self._masterfile(Fin)
    # TODO: consider using pickle 
        #### Requires review should be removed
        align = super().readPhylip()
        self._seqMatchIter(align, master_name, Fout, loc, jsonflag, pangotime)
        if int(self.noiseflag):
            print('Sequence deduplicated. Output is %s\nfor locus %s' % (Fout.name, loc))
#        print('Sequence deduplicated. Output is %s.\nThe number of unique sequenes is %d for locus %s' % (Fout.name, len(idcodedict), loc))
#         super().writejson(idcodedict)
        return

class Load (Taxaswitch):
    #  The noiseflag has been passed through locFasta for phylip files
    #  This could be clarified but it will work bcause of the 1/0 assessment
    #  Startfile is the modified original file post Perl pie
    def __init__(self, locPhylip, locFasta, pathDir, IDdict, phyalignT, pangoFile):
        # TODO: This needs changing and all 'any' need tracking
        if str(locFasta) == '1' or str(locFasta) == '0':
            noiseflag = locFasta
            super().__init__(locPhylip, noiseflag, pathDir, IDdict, phyalignT, pangoFile)
        # if int(locPhylip) == any([0,1]):
        #     self.noiseflag = locPhylip
        self.locFasta = locFasta
        self.pathDir = pathDir
        self.IDdict = IDdict
        self.phyalignT = phyalignT
        self.pangoFile = pangoFile

    def _writephylip(cls, pathDirPhylip, faalign, check=0, locFasta = '', noiseflag=1, rep=1):
        # phylipID_width = cls._checkIDwidth(pathDirPhylip,locFasta)
        phylipID_width = 150
        # What this does is prepare data for parallelisation
        # The distance between the name and the sequence is important in relaxed phylip
        # This approach is dynamic to minimise the distance between the end of the sequence id and start of the sequence
        # The alternative approaches are to just say the largest distance in the input file is the fixed distance
        # For everything - that is too slow. If the number of characters in the input file is not measured it will
        # Trigger bugs in relaxed phylip format 
        phylipList=[]
        # srep = re.search(r'([0-9]+)\.[^.]*\.?[0-9]{0,3}\.fa[a-z]{0,3}',str(locFasta))
        # repNew = str(int(srep.group(1))+1)
        repNew = str(int(rep)+1)
        # phyfile = re.sub(r'([0-9]+)(\.[^.]*\.?[0-9]{0,3}\.)fa[a-z]{0,3}', repNew + r'\2' + 'phy', str(locFasta))
        phyfile = re.sub(r'\.?[0-9]{0,3}\.fa[a-z]{0,3}', r'.' + repNew + '.phy', str(locFasta))
        # isinstance(check, int) is also viable
        if check > 0:
            # phywrite = Path(pathDirPhylip, re.sub(r'\.phy','.' + str(check) + '.phy', phyfile))
            phywrite = Path(pathDirPhylip, re.sub(r'\.phy','.' + str(check) + '.phy', phyfile))
            if int(noiseflag):
                print ('Sub-file has been output, %s' % (Path(phywrite).name))
        phylipList = [phywrite, faalign]
        charaCount = []
        for record in faalign:
            character = int(len(re.findall(r'\S', record.id)) + 1)
            charaCount.append(character)
        phylipID_width = max(charaCount)
        phylipList.append(phylipID_width)
        phylipList.append(rep)
        return phylipList

    def _parallelPrep(cls, pathDirPhylip, faalign, chunksize, locFasta, noiseflag, rep):
        from multiprocessing import Pool, cpu_count
        limit = chunksize + 5000
        if len(faalign) > limit:
            truePhylipList=[]
            processes_count = (cpu_count()-1)
            pool = Pool(processes_count)
            flag = 0
            part = 1
            myrange = list(range(0,len(faalign), chunksize))
            if len(faalign) - myrange[-1] < chunksize:
                myrange.append(len(faalign))
            for j in myrange:
                if j == 0:
                    continue
                tmplist = cls._writephylip(pathDirPhylip, faalign[flag:j,:], part, locFasta, noiseflag, rep)
                truePhylipList.append(tmplist)
                flag = j
                part += 1
            pooledout = pool.map_async(parallel, [row for row in truePhylipList],chunksize=3)
            pool.close()
            pool.join()
            pool.terminate()
            pooledout.get()
            if int(noiseflag):
                print ('This file has been split into %s smaller files of %s sequences per file.' % (part, str(chunksize)))
            raise SystemExit()
        else:
            return 1

    def IDcoding (cls, phyalign, IDdict):
        count = 1
        # if the pangotime flag was mis-set it is reset to zero, i.e. no
        # no pangotype and time parameters
        for record in phyalign:
            base36 = base_repr(count, 36)
            # open SARScov2_5UTR_1-0.25mil.aln.fa
            IDdict[base36] = record.id
            record.id = base36
            record.description = base36
            count += 1
        # print (json.dumps(IDdict, indent = 4))
        return IDdict, phyalign

    def _laurawrite(cls, laurapath, output):
        if laurapath.is_file():
            with laurapath.open("a") as f:
                f.write(output)
            if re.search('failed', output):
                print ("This is the wrong file. The original"
                       " is needed.")
                raise SystemExit('Discontinued.')
        else:
            raise SystemExit("The small logfile called .laura is not available."
                           "Please check %s" % (laurapath.parent))

    def _grep(cls, pattern, rootfile):
        c = re.compile (pattern)
        with io.open(rootfile, "r", encoding="utf-8") as fin:
            if c.search(fin.read()):
                return True
            else:
                return False

    def _perlpie(cls,rootfile):
            print ("Lineage type and/or isolate date analysis is disabled")
            subprocess.call(["perl", "-p", "-i", "-e", "if(/^>.*\|/){s/[ |_]//g}", rootfile], text=True)
            subprocess.call(["perl", "-p", "-i", "-e", "if(/^>.*\|/){s/\|/_/g}", rootfile], text=True)
            subprocess.call(["perl", "-p", "-i", "-e", "if(/^>.*:/){s/:/=/g}", rootfile], text=True)
            subprocess.call(["perl", "-p", "-i", "-e", "if(/^>/){s/[_]{2,3}/_/g}", rootfile], text=True)

    def perlpie (cls, rootfile, pangotime, rep):
        if not str(pangotime) == '1' and not str(pangotime) == '2':
            pangotime = 0
            if not cls._grep:
                cls._perlpie(rootfile)
        if str(rep) == '0' and str(pangotime) == '0':
            subprocess.call(["perl", "-p", "-i", "-e", "s/^(>[^_]+?)_.*/$1/g", rootfile], text=True)
            subprocess.call(["perl", "-p", "-i", "-e", "if(/^>/){s/\.[0-9]+.*$//g}", rootfile], text=True)
            return 1
        laura = Path(rootfile.parent, '.laura')
        if str(pangotime) == '1' and laura.is_file():
            with open(rootfile, 'r') as din:
                for line in din:
                    if line[0] == '>':
                        if re.match(r'^>[^_]+?$', line):
                            cls._laurawrite(laura, "Clean failed pangotime 1")
                        elif re.match(r'^>[^_]+?_[^_]+?$', line):
                            cls._laurawrite(laura, "Already cleaned pangotime 1")
                            break
                        else:
                            cls._laurawrite(laura, "Cleaned\n")
                            break
                    else:
                        print ('The first line of %s is not a fasta header' % (rootfile))
                        raise SystemExit('Discontinued.')
            subprocess.call(["perl", "-p", "-i", "-e", "s/(^>[^_]+)\.[0-9]+(_[^_]+?)_.*/$1$2/g", rootfile], text=True)
            return 1
        if str(pangotime) == '2' and laura.is_file():
            with open(rootfile, 'r') as din:
                for line in din:
                    if line[0] == '>':
                        if re.match(r'^>[^_]+?$', line):
                            cls._laurawrite(laura, "Clean failed pangotime 2")
                        elif re.match(r'^>[^_]+?_[^_]+?$', line):
                            cls._laurawrite(laura, "Clean failed pangotime 2")
                        elif re.match(r'^>[^_]+?_[^_]+?_[^_]+?_[^_]+?_[^_]+?$', line):
                            cls._laurawrite(laura, "Already cleaned pangotime 2")
                            break
                        else:
                            cls._laurawrite(laura, "Cleaned")
                            break
                    else:
                        print ('The first line of %s is not a fasta header' % (rootfile))
                        SystemExit('Discontinued.')
            subprocess.call(["perl", "-p", "-i", "-e", "s/(^>[^_]+)\.[0-9]+(_[^_]+?_[^_]+?_[^_]+?_+[^_]+?)_.*/$1$2/g", rootfile], text=True)
            return 1
    # def getrep (cls, pathDir):
    #     _, rep = re.split('phylip',str(pathDir))
    #     rep = int(rep)
    #     return rep
    

    @classmethod
    def files(cls, locPhylip='full_5UTR_480SARS-CoV-2.phy', locFasta='full_5UTR_480SARS-CoV-2_32_dedup.fa', pathDir='', filechunk= 10000, rep=0, pangotime=0, noiseflag=1, dealDict={}):
        if re.search (r'\.phy$', str(locPhylip)):
            startfile = locFasta
            pangoTmp = re.sub(r'\.fa$', '.v' + pangotime + '.pango.fa', str(startfile))
            parentPath = Path(pathDir).parent
            pangoFile = Path(parentPath, pangoTmp)
            startfile = Path(parentPath,startfile)
            rootfile = Path(pathDir, locPhylip)
            IDdict,phyalignT = {},{}
            print (noiseflag)
            cl = cls(locPhylip, noiseflag, pathDir, IDdict, phyalignT, pangoFile)
            print ('and here too')
        # TODO: switch to pickle
            phyalign = AlignIO.read(rootfile, "phylip-relaxed")
            IDdict, phyalignT = cl.IDcoding(phyalign, IDdict)
            if str(pangotime) == '1' or str(pangotime) == '2':
                if not pangoFile.is_file():
                    subprocess.call(["cp", startfile, pangoFile], text=True)
                    cl.perlpie(pangoFile, pangotime, rep)
            del phyalign
            return cls(locPhylip, noiseflag, pathDir, IDdict, phyalignT, pangoFile)
        # This is leveraging the classmethod to prepare the inherited variables
        # In particular this ensures that only a phylip format is inherited
        # The input is generally raw fasta format and needs the headers preprocessing
        # It is much quicker to do this using perl pie than opening and rewriting
        # The processing time needs to be as quick as possible.
        elif re.search (r'\.fa[a-z]{0,3}$', locFasta):
            pangoFile = ''
            IDdict = {'a':'b'}
            phyalignT = {'c':'d'}
            # cl = cls(locFasta, '', pathDir, IDdict, {})
            cl = cls('', locFasta, pathDir, IDdict, {}, pangoFile)
            # rep = cl.getrep(pathDir)
            if not pathDir.is_dir():
                os.mkdir(pathDir)
            # The two directory approach is explained through the code
            rootfile = Path(Path(pathDir).parent, locFasta)            
            # Perl pie is excuted on the parent directory hence 'fawrite'
            # Perl pie was faster than using this within python,
            if str(rep) == '0':
                locFasta_copy = re.sub(r'\.fa$', '.taxaid.fa', locFasta)
                rootfile_copy = Path(Path(pathDir).parent, locFasta_copy)
                # copy2(rootfile, rootfile_copy)
                subprocess.call(["cp", rootfile, rootfile_copy], text=True)                
                # cl.perlpie(rootfile, pangotime, rep)
                cl.perlpie(rootfile_copy, pangotime, rep)
            # Switch to read pickle
                faalign = AlignIO.read(rootfile_copy, "fasta")
            else:
            # Switch to read pickle
                faalign = AlignIO.read(rootfile, "fasta")
            cl._parallelPrep(pathDir, faalign, filechunk, locFasta, noiseflag, rep)
            # The checklength routine isn't well named but gives a gist.
            # If the file size is bigger than 5000 + optimised chunk 1750
            # checklength will terminate following splitting the file, otherwise it continues
            # this time the output file will not be dumped in the 'phylip'rep no. directory
            # but in the  parent directory
            locPhylip = cl._writephylip(pathDir, faalign, 1)
            # The code above is needed to stop a file of the required size being rewritten
            locFasta = 'na'
            return cls(locPhylip.name, locFasta, pathDir, IDdict, phyalignT, '')

    @staticmethod
    def checkit(phylipFile):
        return AlignIO.read(phylipFile, "phylip-relaxed")

    @staticmethod
    def debugit(pathdir, fastafile):
        results = []
        cnt=-1
        mypath = Path(pathdir.parent, fastafile)
        with open(mypath, 'rt') as t:
            line = t.readlines()
            length = len(line) - 1
            while line:
                if cnt == length:
                    break
                cnt+=1
                if line[cnt][0] == '>':
                    if len(line[cnt]) > 50:
                        results.append(len(re.findall(r'\S', line[cnt]))+1)
        if len(results) > 0:
            print (f"{'Largest taxa label is %s' % (str(max(results)))}")
        else:
            results.append(50)
        return max(results)

    @staticmethod
    def bigGenotype (outfile):
        fout = subprocess.call(["perl", "-n", "-l", "-e", "print $1 if (m/([A-Za-z0-9\|]+)\s[A-Z\-]/g);", outfile], text=True)
        return fout

# Writing out the Phylip files is parallelised and this occurs in global space.
# There is no other way to due this in parallelisation. If it is single threaded
# then it will take place with the classmethod
def parallel(inputList):
    phywrite = inputList[0]
    faalign = inputList[1]
    phylipID_width = inputList[2]
    rep = inputList[3]
    if str(rep) == '0':
        faalign = delNs(phywrite, faalign) 
    with open(phywrite, 'w') as output_handle: # phylipID_width
        SequentialPhylipWriter(output_handle).write_alignment(faalign, id_width=phylipID_width)
    # TODO:  Add here dump to pickle

def delNs (phywrite, faalign):
    locFasta = re.sub(r'\.[^\.]+\.[^\.]+\.[^\.]+\.[^\.]+$','',phywrite.name)
    pathDir = phywrite.parent.parent
    lotsNs = re.compile('[N]{20}')
    alignWithoutNs = Bio.Align.MultipleSeqAlignment([])
    removedIDs = []
    for record in faalign:
        if lotsNs.search(str(record.seq)):
            removedIDs.append(record.id)
            continue
        else:
            alignWithoutNs.append(record)
    info = delNstring(removedIDs)
    del faalign
    thedate = '\n' + datetime.today().strftime('%Y-%m-%d %H:%M:%S')
#        print (info + thedate)
    # m = re.match(r'(.+)\.+',locFasta)
    # rmfout=  m.group(1) + '.removedSeq.txt'
    rmfout=  locFasta + '.removedSeq.txt'
    with open(Path(pathDir, rmfout), "a") as f:
        f.write(str(info) + thedate)
    return alignWithoutNs

def delNstring (removedIDs):
    return f"{'Removed %s taxa due to >24 continuous Ns, viz. %s' % (len(removedIDs), ' '.join([a for a in removedIDs]))}"
