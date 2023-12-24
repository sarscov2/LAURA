#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 07:43:47 2022
@author: michaelgaunt
""" 

# from trickstest import Convert, Input

from biolaura2 import Biotricks as Bt
from dblaura import ParellelFAio as Pfaio
from rnalaura2 import Rnatricks as Rnat, Rnaloop, Alignstruc
from simlaura2 import Shufflesim
from qclaura2 import RnaQC
from yesnolaura import Yesno
from taxatricks24d_8g_2 import Load
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from pathlib import Path
import json, re, os, time, shutil, gc, copy, argparse, errno, sys #, keyboard
from collections import Counter, defaultdict
from random import seed
from Bio import AlignIO

processes_count = (cpu_count()-1)
# pathwayTrue goes directly to the phylip directory
pathwayTrue = ''
deduplicat = 'dedup.fa'
allelecat = 'allele.fa'
deduplicatJson = 'dedup.json'

def parallize (inputFile):
    # NOTE the input here is inputFile.name ... its pathlib f
    # so the resulting string is everything after the final /
    idDB = inputFile[1]
    loci, jsonflag, alleleflag, noiseflag, pangotime, startfile, hdf5flag, rnamode = re.split (r'[\|]', str(inputFile[0].name))
    originalFile = Path(inputFile[0].parent.parent.parent, startfile)
    myfiles4 = Load.files(str(Path(inputFile[0].parent).name), originalFile, str(inputFile[0].parent.parent), 1750, 1, pangotime, noiseflag, hdf5flag, {}, idDB, rnamode, 0)
    # Zero refers to zerostartmode which isn't relevant
    return myfiles4.seqMatch(loci, jsonflag, alleleflag, pangotime, hdf5flag)

def startit(dealDict, size, rep, pangotime, noiseflag, hdf5flag, fastaLoc = '.fa', rnamode = 0, zerostartmode=0):
# Note that size is 'filechunk' within taxatricks
    # idDB = {}
    return Load.files('', fastaLoc, pathwayTrue, size, rep, pangotime, noiseflag, hdf5flag, dealDict, {}, rnamode, zerostartmode)

def rnaloopsim(myseed, pathway, rnaref):
    records = SeqIO.parse(Path(pathway,rnaref),'fasta')
    iniSs = Shufflesim()
    simout = iniSs.getsimout(pathway, rnaref)
    del iniSs
    rnastruct = ''
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
        iniSs.writesim(1)
        iniSs.writestems()
        myseed  += 1 

def genCoM(locusVals,myloops, pathway, rnaref):
    CoM = []
    if re.search('-', str(locusVals)):
        CoM.append(locusVals)
    else:
        loopsearch = re.split(r'[ ,;:]', myloops)
        iniLoop = Rnaloop(pathway, rnaref)
        loci,score = iniLoop.getloops()
        CoM = [loci[int(l) - 1] for l in loopsearch]
        CoM = [str(tup[0]) + '-' + str(tup[1]) for tup in CoM]
    return CoM

def structuralAlign(rnaalignflag, startfile, pathway):
    start = time.time()
    if rnaalignflag:
    # TODO: This needs changing structuralAlign
        junktest(startfile)
        iniAl = Alignstruc(pathway, startfile, 1, 100)
        # iniAl = Alignstruc(pathway, firstfile, 1, 10)
        iniAl.rnaParallel()
        end = time.time()
        print(
            '{} execution took {} seconds.'.format(
                'RNA structure alignment', round(end-start)
                )
            )

def parallizeKitty(kitty):
    dealedDict_tmp = {}
    count = kitty[0]
    for k,v in kitty[1].items():
        n = int(((count) %20))
        newFout = str(n)
        if k:
            if newFout in dealedDict_tmp:
                dealedDict_tmp[newFout].append([k, v])
            elif not newFout in dealedDict_tmp:
                dealedDict_tmp[newFout] = [[k, v]]
    return dealedDict_tmp

def alleleKitty(kitty):
    return AlignIO.read(kitty[1], 'fasta')

def filechecking(rootname, tracker, info, hdf5flag,rnamode):
    thisrep = str(info['rep'])
    thisrepLessOne = str(info['rep'] - 1)

    if pathwayTrue.is_dir():
        pass
    else:
        if str(rnamode) == '1':
            setpathway(str(pathwayTrue.parent), 'rna' + thisrepLessOne)
        else:
            setpathway(str(pathwayTrue.parent), 'phylip' + thisrepLessOne)
    if hdf5flag: 
        rootlist = pathwayTrue.glob(rootname + '*' + '.hdf5')
    else:
        rootlist = pathwayTrue.glob(rootname + '*' + '.phy')
    rootfileBool = []
    
    for r in rootlist:
        rootfileBool.append(Path(r).is_file())
    rootfileBool = all(rootfileBool)
    return rootfileBool,thisrep,thisrepLessOne

def getfiles(startfile, pathway):
    filesuffix= re.compile(r'[0-9]{0,3}.fa[a-z]{0,3}')
    if filesuffix.search(startfile):
        rootname, _ = filesuffix.split(startfile)
    else: 
        rootname, _ = re.split(r'.fa[a-z]{0,3}$',startfile)
    # Sometimes the file 0.fa will be present if there is no associated .rep file 
    # it will then redo the phylip breakup and there is no need to rename the file
    oldlist = Path(pathway).glob('*' + '.' + '0' + '.fa')
    oldfile = [o for o in oldlist]
    oldfile.append(Path('junk'))
    oldfile = oldfile.pop(0)
    return oldfile, rootname

def loadparallel (rootname, locus, jsonflag, alleleflag, noiseflag, pangotime, startfile, hdf5flag, rnamode, full):
# TODO: Bug
    if hdf5flag:
        filePath = pathwayTrue.glob(rootname + '*.?*.hdf5')
    else:
        if str(rnamode) == '1' and full:
            if len(Path(rootname).parents) > 0:
                rootname = str(Path(rootname).name)
        filePath = pathwayTrue.glob(rootname + '*.?*.phy')
    varload = locus + '|' + str(jsonflag) + '|' + str(alleleflag) + '|' + str(noiseflag)  + '|' + str(pangotime) + '|' + str(startfile) + '|' + str(hdf5flag)  + '|' + str(rnamode)
# # The async_map is converted into a multivariable function by using Path
    files = [Path(f, varload) for f in filePath]
    return files

def getdeduplicated(suffix, allelesuffix, info, pathway, rootname, alleleflag,locus, pangotime):
    alleleScoop = []
    foutAllele = ''
    # filescoop = pathwayTrue.glob(Path(rootname).name + '*' + suffix)
    # if re.search('dedup', rootname):
        # rootname = re.sub(r'(\.?[0-9]*){0,2}\.dedup.fa', '',rootname )
    Fout = Path(pathway, rootname + '.' + str(info['rep']+1) + '.' + suffix)
    if str(alleleflag) == '1' or str(alleleflag) == '2':
        alleleScoop = pathwayTrue.glob(Path(rootname).name + '*' + allelesuffix)
        rootname = re.sub(r'\.dedup.fa', 'locus' + locus + '.allele.fa', rootname )
        foutAllele = Path(
            pathway, Path(rootname)
            .name + '.' + 'v' + str(pangotime) + '.' + 'rep' +  str(info['rep']+1) + '.' + locus + '.' + allelesuffix
            )
        print ('Convergence is imminent.\nThe raw output is, %s' % (
            Path(Fout).name))
    return Fout, alleleScoop, foutAllele

def fastadumpParallel(seqList, allele=0):
    dealDict = defaultdict(list)
    seqList2 = [[count, seqkv] for count, seqkv  in enumerate(seqList)]
    del seqList
    pool2 = Pool(processes_count)
    if str(allele) == '0':
        poolout2 = pool2.map_async(parallizeKitty, [kitty for kitty in seqList2], chunksize=5)
        iniRna = Rnat('','',pool2)
        pooloutGet = iniRna.pooladmin(poolout2) # check this parameter
        for wrecord in pooloutGet:
            for li in wrecord[''.join(wrecord.keys())]:
                dealDict[''.join(wrecord.keys())].append(li)
    elif str(allele) == '1' or str(allele) == '2':
        poolout2 = pool2.map_async(alleleKitty, [kitty for kitty in seqList2], chunksize=2)
        iniRna = Rnat('','',pool2)
        pooloutGet = iniRna.pooladmin(poolout2) 
        alignList = [record for records in pooloutGet for record in records]
        iniBt = Bt()
        dealDict = iniBt.makealign(alignList)
    else:
        return False
    del pooloutGet, poolout2, seqList2, iniRna
    gc.collect()
    return dealDict

def jsondump(catDedupJson, FoutJson):
    resultjson = []
    for kitty2 in catDedupJson:
        with open(kitty2, "r") as djin:
            # print (kitty2)
            resultjson.append(json.load(djin))
    with open(FoutJson, "w") as djout:
        json.dump(resultjson, djout)

def checkphylip (pathway = '', folder=''):
    startfolder = Path(pathway, folder)
    ch = startfolder.glob('*' + '.phy')
    check = [c for c in ch] 
    if len(check) > 1:
        tracker.touch()
    else:
        raise Exception ('I cannot find the base phylip files')
        return False

def startlaura(startfile, pathway,tracker, size, pangotime,noiseflag,hdf5flag,rnamode, zerostartmode, tab):
    # oldfile here is a 0.fa file then pick it up and use it
    startPath = Path(pathway,startfile)
    if startPath.is_file() == True and not Path(pathway, '.laura').is_file():
        try:
            startit({}, size, 0, pangotime, noiseflag, hdf5flag, startPath, rnamode, zerostartmode)
        except SystemExit:
            if str(rnamode) == '1':
                startdir = 'rna0'
            else:
                startdir = 'phylip0'
            checkphylip (pathway,startdir)
            return True
    # tab is in operation
    elif str(tab) == '2':
        return True
    else:
        print (
            'I cannot find either the directory {} or the starting file.\n{}'
            .format(pathwayTrue.name,'If a complete restart is needed please delete file .laura')
            )
        return False

def junktest(firstfile):
    if 'junk' in str(firstfile):
        raise SystemExit(
            'Cannot find the file ending with {}.\nRun nucleotide alignment first'
            .format('0.fa'))

def preplaura(pathway, info):
    laurapath = Path(pathway, '.laura')
    if laurapath.is_file() and info['rep'] == 0:
        laurapath.unlink()

def makefirstfile(rnamode, info, startfile, pathway):
    if str(rnamode) == '1':
        mydir = 'rna' + str(info['rep'])
        firstfiletmp, rootname =  getfiles(startfile, pathway)
        junktest(firstfiletmp)
        firstfile = Path(re.sub('fa$','rnaStr.fa',str(firstfiletmp)))
    else:
        mydir = 'phylip' + str(info['rep'])
        firstfile, rootname =  getfiles(startfile, pathway)
    return firstfile, rootname, mydir

def inforep(rnamode, info, pathway):
    repCount = False
    if str(rnamode) == '1':
        mydir = 'rna' + str(info['rep'])
        repCount = True        
    else:
        mydir = ['phylip' + str(int(info['rep']) - replicate) for replicate in range(-1,3)]
        if Path(pathway,mydir[1]).exists():
            info['rep'] += 1
            mydir = mydir[0]
        elif Path(pathway,mydir[2]).exists():
            mydir = mydir[1]
            repCount = True
        elif Path(pathway,mydir[3]).exists():
            raise RuntimeError(
                'The directory {} is out of step with replicate no.'.format(
                    Path(pathwayTrue,mydir[2])))
        else:
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 
                str(mydir[0] + ' or ' + str(mydir[1]))
                )
    return mydir, repCount

def argparsing(startfile, pathway):
    parser = argparse.ArgumentParser(
        description='Assess haplotype diversity for a single RNA secondary structure',
        epilog='qc = quality control; sim = simulation rnastruc; rnaalign = RNA structural alignment'
        )
    parser.add_argument(
        '--input', action='store', default=startfile, required=False, 
        help="the input file must be fasta format")
    parser.add_argument(
        '--qc', action='store_true', default=False, required=False, 
        help="perform a QC on the input fasta file")
    parser.add_argument(
        '--sim', action='store_true', default=False, required=False, 
        help="perform a simulation using a reference file")
    parser.add_argument(
        '--rnaalign', action='store_true', default=False, required=False, 
        help="produce a RNA alignment from a fasta file")
    args = parser.parse_args()
# TODO: to change later
    if args.input and not pathway:
        if str(Path(args.input).parent) == '.':
            pathway = Path(__file__).absolute().parent
        else:
            pathway = Path(args.input).parent
    file = Path(args.input).name
    if not Path(pathway, file).exists():
        print ("Unable to find {}".format(Path(pathway,file)))
        sys.exit(0)
    return file, pathway, args

def setpangoAllele(alleleflagset, pangotimeflag):
    if info['terminate'] > 0.7 or info['terminate'] == 0.49:
        alleleflag = alleleflagset
        pangotime = pangotimeflag
        return pangotime, alleleflag
    else:
        return 0, 0

def rnastrucMode (startfile, pathway):
    pcheck = Path(pathway, startfile)
    # try:
    if pcheck.exists():
        rnafile = Path(re.sub('0.fa','rnaStr.fa', str(pcheck)))
        if rnafile.exists():
            return rnafile
        else:
            altfile = list(pcheck.parent.glob('*' + '.rnaStr.fa'))            
            if len(altfile) == 1:
                question = "I've found a file called {}.\n Is this the file for\
                      RNA structure analysis?".format(altfile[0])
                iniYesno = Yesno(question)
                answer = iniYesno.qResponse()
                if answer:
                    return str(altfile[0])
                else:
                    print ("Okay, please enter restart LAURA using \
                           '--structure' or '-s' mode naming the structure\
                               file.")
                    sys.exit(0)
            else:
                print("The file containing the RNA structure data cannot be \
                      found.\n It should end with 'rnaStruc.fa'. The file \
                          could be absent, or there could be muliple \such files")
                sys.exit(0)
    # except:
        # raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), 
                                # str(startfile))

def setpathway(pathway, mydir):
    global pathwayTrue
    pathwayTrue = Path(pathway,mydir)

if __name__ == '__main__':

    startfile, rnaref = 'SARScov2_5UTR_1-0.75mil.full.aln2.qc.0.fa', 'rnaref.fa'
    pathway = "/Volumes/Data/5UTR_millcuarto/fullTest/fullTest2/"
    file, pathway, args = argparsing(startfile, pathway)
    locusVals, myloops = '1-400', '1'
    full = True
    info = Counter()
    info, lastdict = {'terminate': 0, 'rep': 0, 'restart':0}, {}
    replicants = 200
    rootfile,oldfile = '', ''
    hdf5flag, jsonflag = 0, 0
    alleleflag, alleleflagset = 0, 1    # needs automatically assigning
    noiseflag, tab = 0, 0
    pangotime, pangotimeflag = 0, 0 # this will be 0 (nothing), 1 (pangotype)
                                    # controls the pangotime within the last loop. 
                                    # Flags better combined. 1 or 2 (pangotype and time)
    writepango = 0 #  writes out to disk during parallelisation
    rnastruc, rnamode = 0, 0 # high level control. This is a static variable identity     
    zerostartmode = 1 # If the 0.fa files is generated this will start from that point
    size, chunk = 1750, 1

    if full:
        rnafile = rnastrucMode(startfile, pathway)
        rnastruc = 1

    args.rnaalign = False    

    tracker = Path(pathway, '.laura')
    if args.sim:
        myseed = 10
        rnaloopsim(myseed, pathway, rnaref)
    elif args.qc:
        print('QC commencing ...') 
        qcfile = Path(pathway,startfile)
        iniQC = RnaQC(qcfile)
        fout, _, _, _ = iniQC.coreQC(5, locusVals)
    elif args.rnaalign:
        structuralAlign(args.rnaalign, startfile, pathway)
    else:
        CoM = genCoM(locusVals,myloops, pathway, rnaref)
    if not args.sim and not args.qc and not args.rnaalign:
        # Section 0: loc CoM
        if full:
            rnaterminate = 3
            firstreprna = True
        else:
            rnaterminate = 1
        while rnastruc < rnaterminate:
            print (rnastruc, "rnastruc variable")
            if rnastruc > 1 and firstreprna:
                startfile = rnafile
                info['rep'] = 0
                info['terminate'] = 0
                rnamode = 1
                firstreprna = False
                if Path(tracker).exists():
                    os.unlink(tracker)
            for locus in CoM:
                # Section 0, full flag
                preplaura(pathway, info)
                while info['terminate'] < 1 and info['rep'] < replicants:
                    start = time.time()
                    firstfile, rootname, mydir = makefirstfile(
                        rnamode, info, startfile, pathway
                        )
                    # TODO: Bug rootname
                    setpathway(pathway, mydir)
                    rootfile = Path(
                        pathway,rootname + '.' + str(info['rep']) + '.fa'
                        )
                    # got remove .laura if the replication is out of step with 
                    # master file remove all other files 
                    print (info['rep'])
                    if tracker.is_file() and info['rep'] > 0:
                        rootfileBool,thisrep,thisrepLessOne = filechecking(
                            rootname, tracker, info, hdf5flag, rnamode
                            )
                        if not rootfileBool:
                            print (
                                'Cannot to find the input phylip files. Terminated.'
                                )
                            break
                    # Both the tracker .rep must be present and the next base fasta file it specifies
                    elif info['rep'] == 0:
                        if str(rnamode) == '1':
                            tab += 1
                            startfile = firstfile.name
                            zerostartmode = 0
                        checkcont = startlaura(
                            startfile, pathway, tracker, size, pangotime,
                            noiseflag,hdf5flag,rnamode,zerostartmode,tab
                            )
                        if checkcont is False:
                            break
                    else:
                        print ("Can not find the start file or restart file")
                        break
                    # SECTION 1
                    try:
                        # print(alleleflag, 'alleleflag')
                        files = loadparallel(
                            rootname, locus, jsonflag,alleleflag, noiseflag, 
                            pangotime, startfile, hdf5flag, rnamode, full
                            )
                        if len(files) == 0:
                            raise Exception()
                    except:
                        raise Exception (
                            "I can not file the phylip files in %s" % (mydir)
                            )
                    # # SECTION 1
                    idDB = {}
                    if pangotime > 0:
                        print ('Genetic and epi info DB initiating')
                        iniPfaio = Pfaio(
                            Path(
                                pathway,startfile),pangotime, info['rep'], args.qc, writepango
                            )
                        idDB = iniPfaio.parallDictWrit()
                    pool = Pool(processes_count)
                    poolout = pool.map_async(
                        parallize,[[row, idDB] for row in files], chunksize=chunk
                        )
                    iniPool = Rnat('','',pool)
                    poolOutget = iniPool.pooladmin(poolout)
                    seqList = []
                    for seqkv in poolOutget:
                        for k,v in seqkv.items():
                            seqList.append({k:v})
                        # seqDict[k] = tuple(seqDict[k] for seqDict in poolOut.get()) 
                    del iniPool
                    # # SECTION 2
                    pangotime, alleleflag = setpangoAllele(
                        alleleflagset, pangotimeflag
                        )
                    foutDedup, alleleDedup, foutAllele = getdeduplicated(
                        deduplicat, allelecat, info, pathway, rootname, alleleflag, locus, pangotime
                        )

                    # a +1 has been added to this file 
                    # this is formalised just before the convergence check
                    dealDict = fastadumpParallel (seqList)
                    
                    if str(alleleflag) == '1' or str(alleleflag) == '2':
                        alleleDict = fastadumpParallel (
                            alleleDedup, alleleflag
                            )
                    if jsonflag == 1:
                        catDedupJson, foutDedupJson = getdeduplicated(
                            deduplicatJson,info,'', pathway, rootname, ''
                            )
                        jsondump (catDedupJson, foutDedupJson)
                    # SECTION 3
                    # Have to reset the path so the new phylip files go into a new directory
                    # +1 is essential to all the info['rep] to distinguish it 
                
                    mydir, repCount = inforep(rnamode, info,pathway)
                    
                    setpathway(pathway, mydir)
                    iniBt = Bt(dealDict)
                    currentSeqDict = iniBt.reformShuffle()
                    if info['rep'] > 2:
                        print (
                            str(len(lastdict)) + ' ' + str(
                                len(currentSeqDict))
                            )
                        if len(lastdict) == len(currentSeqDict):
                            if lastdict == currentSeqDict:
                                # info['terminate'] += 0.4
                                info['terminate'] += 0.1
                        lastdict = copy.copy(currentSeqDict)
                    bioalign = iniBt.bioobj()            
                    # print ( 'terminate score', info['terminate'])
                    try:
                        if info['terminate'] < 1:
                            zerostartmode = 0
                            startit(
                                bioalign, size, info['rep'], pangotime, noiseflag, 
                                hdf5flag, firstfile, rnamode, zerostartmode
                                )
                        else:
                            alignFin = iniBt.makealign(bioalign)
                            if re.search('\.\.', str(foutDedup)):
                                foutDedup = Path(re.sub('\.\.', '.', str(foutDedup)))
                                foutAllele = Path(re.sub('\.\.', '.', str(foutAllele)))
                            AlignIO.write(alignFin, Path(pathway,foutDedup), 'fasta')
                            AlignIO.write(alleleDict, Path(pathway,foutAllele), 'fasta')
                            break
                    except SystemExit:
                        del dealDict, currentSeqDict, iniBt, bioalign
                        gc.collect()

                        if len(files) == 1:
                            info['terminate'] += 0.49
                            print (
                                "Less than 1750 sequences remain in deduplication\nReplication no. is %s\n" % (
                                    info['rep'])
                                )
                        print ('Replicant %s completed' % (info['rep']))
                        if repCount:
                            info['rep'] += 1
                        if str(rnamode) == '1':     
                            dirOld = Path(pathway, 'rna' + str(info['rep']-2))
                            if dirOld.is_dir():
                                shutil.rmtree(dirOld)
                        dirOld = Path(pathway, 'phylip' + str(info['rep']-2))
                        if dirOld.is_dir():
                            shutil.rmtree(dirOld)
                        # Needs to be a list because might be empty
                        continue
                if info['terminate'] < 1:
                    print(
                        "LAURA has not converged but %s replications are complete" 
                        % (info['rep']))
                print ( 'terminate score', info['terminate'])
                if info['terminate'] >= 1:
                    print ("LAURA has converged at replicate %s." % (info['rep'])) 
                    if str(rnamode) == '1':
                        dirOld = Path(pathway, 'rna' + str(info['rep']-2))
                    else:
                        dirOld = Path(pathway, 'phylip' + str(info['rep']-2))
                    if dirOld.is_dir():
                        shutil.rmtree(dirOld)
                end = time.time()
                print(
                    '{} execution took {} seconds.'.format('Total', round(end-start)))
            rnastruc += 1
        # print ('RNA structural alignment identity is now complete')
    else:
        if args.qc:
            print ("QC is complete")
        elif args.sim:
            print ("Simulation file output")