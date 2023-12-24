#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 19:57:10 2022

@author: michaelgaunt

Code integrated into the run file
"""

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

def fastadumpParallel(seqList, allele=0):
    dealDict = defaultdict(list)
    seqList2 = [[count, seqkv] for count, seqkv  in enumerate(seqList)]
    del seqList
    pool2 = Pool(processes_count)
    if allele == 0:
        poolout2 = pool2.map_async(parallizeKitty, [kitty for kitty in seqList2], chunksize=5)
        iniRna = Rnat('','',pool2)
        pooloutGet = iniRna.pooladmin(poolout2) # check this parameter
        for wrecord in pooloutGet:
            for li in wrecord[''.join(wrecord.keys())]:
                dealDict[''.join(wrecord.keys())].append(li)
    elif allele == 1 or allele == 2:
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


    # return Fout
        # This code was time tested and found to be 3 seconds slower
        #     resultdedup.append(din.readlines())
        # with open (Fout, "w") as dout:
        #     for line in resultdedup:
        #         for l in line:
        #             dout.write('%s' % (l))