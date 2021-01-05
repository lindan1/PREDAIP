# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 09:45:22 2020

@author: 20170426-3
"""
import numpy as np
import math
from collections import Counter


#AAC
def AAC(sequence):
    dataset = []
    for pep in sequence:
        amino = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}
        for i in amino.keys():
            for j in pep:
                if j == i:
                    amino[i] += 1
                else:
                    continue
            amino[i] = amino[i] / len(pep)     
        dataset.append(list(amino.values())) 
    return dataset


#DPC
def DPC(sequence):
    dataset = []
    amino = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    b = []
    for i in amino:
        for j in amino:
            b.append(i + j) 
    for pep in sequence:
        DPC = [0] * 400
        for m, n in zip(pep[:-1], pep[1:]):
            k = m + n
            index = b.index(k)
            DPC[index] += 1 
        dataset.append(DPC)
    return dataset


#CKSAAP
def CKSAAP(sequence, gap=4):

	AA = 'ACDEFGHIKLMNPQRSTVWY'
	dataset = []
	aaPairs = []
	for aa1 in AA:
		for aa2 in AA:
			aaPairs.append(aa1 + aa2)

	for seq in sequence:
		code = []
		for g in range(gap+1):
			myDict = {}
			for pair in aaPairs:
				myDict[pair] = 0
			sum = 0
			for index1 in range(len(seq)):
				index2 = index1 + g + 1
				if index1 < len(seq) and index2 < len(seq) and seq[index1] in AA and seq[index2] in AA:
					myDict[seq[index1] + seq[index2]] = myDict[seq[index1] + seq[index2]] + 1
					sum = sum + 1
			for pair in aaPairs:
				code.append(myDict[pair] / sum)
		dataset.append(code)
	return dataset


#TAAI
import sys, os, re, platform
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)

def TAAI(sequence, use_pro = [87, 379, 29, 88, 145, 30, 25, 26, 27, 95, 452, 92, 154, 15, 81]): 
    fileAAindex = re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + r'codes\AAindex.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + 'codes/AAindex.txt'

    f = open(fileAAindex)
    AAindex_file = f.readlines()  
    AAindex = []
    for i in AAindex_file:
        AAindex.append(i.rstrip().split()[1:] if i.rstrip( ) != '' else None)
    
    amino = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X']
        
    
    dataset = []
    for seq in sequence:
        code = []
        for pro in use_pro: 
            propertyi = AAindex[pro]   
            data = []
            mean_var_std = []
            for aa in seq:
                index = amino.index(aa)
                data.append(float(propertyi[index]))
            data = np.array(data)
            data_mean = np.mean(data)
            data_var = np.var(data)
            data_std = np.std(data)
            mean_var_std.extend([data_mean, data_var, data_std])
            code = code + mean_var_std
          
        dataset.append(code)
        
    return dataset


# DDE
def DDE(sequence):
	AA = 'ACDEFGHIKLMNPQRSTVWY'

	myCodons = {'A': 4, 'C': 2, 'D': 2, 'E': 2, 'F': 2, 'G': 4, 'H': 2, 'I': 3, 'K': 2, 'L': 6, 'M': 1, 'N': 2,
		         'P': 4, 'Q': 2, 'R': 6, 'S': 6, 'T': 4, 'V': 4, 'W': 1, 'Y': 2}

	dataset = []
	diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]

	myTM = []
	for pair in diPeptides:
		myTM.append((myCodons[pair[0]] / 61) * (myCodons[pair[1]] / 61))

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i

	for seq in sequence:
		tmpCode = [0] * 400
		for j in range(len(seq) - 2 + 1):
			tmpCode[AADict[seq[j]] * 20 + AADict[seq[j+1]]] = tmpCode[AADict[seq[j]] * 20 + AADict[seq[j+1]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]

		myTV = []
		for j in range(len(myTM)):
			myTV.append(myTM[j] * (1-myTM[j]) / (len(seq) - 1))

		for j in range(len(tmpCode)):
			tmpCode[j] = (tmpCode[j] - myTM[j]) / math.sqrt(myTV[j])

		dataset.append(tmpCode)
	return dataset




def GAAC(sequence):
    group = {'alphatic': 'GAVLMI',  'aromatic': 'FYW',  'postivecharge': 'KRH',  'negativecharge': 'DE',  'uncharge': 'STCPNQ'}
    groupkey = group.keys()
    
    dataset = []
    for seq in sequence:
        code = []
        count = Counter(seq)
        myDict = {}
        
        for key in groupkey:
            for aa in group[key]:
                myDict[key] = myDict.get(key, 0) + count[aa]
        for key in groupkey:
            code.append(myDict[key]  / len(seq))
        dataset.append(code)
    
    return dataset



#CTDC
def Count1(seq1, seq2):
	sum = 0
	for aa in seq1:
		sum = sum + seq2.count(aa)
	return sum


def CTDC(sequence):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}
	
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	dataset = []
	
	for seq in sequence:

		code = []
		for p in property:
			c1 = Count1(group1[p], seq) / len(seq)
			c2 = Count1(group2[p], seq) / len(seq)
			c3 = 1 - c1 - c2
			code = code + [c1, c2, c3]
		dataset.append(code)
	return dataset




def Count2(aaSet, seq):
	number = 0
	for aa in seq:
		if aa in aaSet:
			number = number + 1
	cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
	cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

	code = []
	for cutoff in cutoffNums:
		myCount = 0
		for i in range(len(seq)):
			if seq[i] in aaSet:
				myCount += 1
				if myCount == cutoff:
					code.append((i + 1) / len(seq) * 100)
					break
		if myCount == 0:
			code.append(0)
	return code


#CTDT
def CTDT(sequence):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}

	
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	dataset = []

	for seq in sequence:
		code = []
		aaPair = [seq[j:j + 2] for j in range(len(seq) - 1)]
		for p in property:
			c1221, c1331, c2332 = 0, 0, 0
			for pair in aaPair:
				if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
					c1221 = c1221 + 1
					continue
				if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
					c1331 = c1331 + 1
					continue
				if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
					c2332 = c2332 + 1
			code = code + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]
		dataset.append(code)
	return dataset


def CTDD(sequence):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}


	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	dataset = []

	for seq in sequence:
		code = []
		for p in property:
			code = code + Count2(group1[p], seq) + Count2(group2[p], seq) + Count2(group3[p], seq)
		dataset.append(code)
	return dataset


#GDPC
def GDPC(sequence):
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	groupKey = group.keys()
	dipeptide = [g1+'.'+g2 for g1 in groupKey for g2 in groupKey]

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	dataset = []

	for seq in sequence:
		code = []
		myDict = {}
		for t in dipeptide:
			myDict[t] = 0

		sum = 0
		for j in range(len(seq) - 2 + 1):
			myDict[index[seq[j]]+'.'+index[seq[j+1]]] = myDict[index[seq[j]]+'.'+index[seq[j+1]]] + 1
			sum = sum +1

		if sum == 0:
			for t in dipeptide:
				code.append(0)
		else:
			for t in dipeptide:
				code.append(myDict[t]/sum)
		dataset.append(code)

	return dataset



#GTPC
def GTPC(sequence):
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	groupKey = group.keys()
	triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	dataset = []

	for seq in sequence:
		code = []
		myDict = {}
		for t in triple:
			myDict[t] = 0

		sum = 0
		for j in range(len(seq) - 3 + 1):
			myDict[index[seq[j]]+'.'+index[seq[j+1]]+'.'+index[seq[j+2]]] = myDict[index[seq[j]]+'.'+index[seq[j+1]]+'.'+index[seq[j+2]]] + 1
			sum = sum +1

		if sum == 0:
			for t in triple:
				code.append(0)
		else:
			for t in triple:
				code.append(myDict[t]/sum)
		dataset.append(code)

	return dataset


#CKSAAGP
def generateGroupPairs(groupKey):
	gPair = {}
	for key1 in groupKey:
		for key2 in groupKey:
			gPair[key1+'.'+key2] = 0
	return gPair


def CKSAAGP(sequence, gap=5):
    
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'}

	AA = 'ARNDCQEGHILKMFPSTWYV'

	groupKey = group.keys()

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	gPairIndex = []
	for key1 in groupKey:
		for key2 in groupKey:
			gPairIndex.append(key1+'.'+key2)

	dataset = []

	for seq in sequence:
		code = []
		for g in range(gap + 1):
			gPair = generateGroupPairs(groupKey) #dictionary
			sum = 0
			for p1 in range(len(seq)):
				p2 = p1 + g + 1
				if p2 < len(seq) and seq[p1] in AA and seq[p2] in AA:
					gPair[index[seq[p1]]+'.'+index[seq[p2]]] = gPair[index[seq[p1]]+'.'+index[seq[p2]]] + 1
					sum = sum + 1

			if sum == 0:
				for gp in gPairIndex:
					code.append(0)
			else:
				for gp in gPairIndex:
					code.append(gPair[gp] / sum)

		dataset.append(code)

	return dataset


 
def ExtractFeatures(sequence, *var): #__file__ is the AAindex filepath
    n = len(var)
    feature = {}
    for i in range(n):
        if var[i] == 1:
            feature[i] = DDE(sequence)
        if var[i] == 2:
            feature[i] = AAC(sequence)
        if var[i] == 3:
            feature[i] = CKSAAGP(sequence, gap=5)
        if var[i] == 4:
            feature[i] = GDPC(sequence)
        if var[i] == 5:
            feature[i] = CTDC(sequence)
        if var[i] == 6:
            feature[i] = CKSAAP(sequence, gap=5)
        if var[i] == 7:
            feature[i] = DPC(sequence)
        if var[i] == 8:
            feature[i] = CTDD(sequence)
        if var[i] == 9:
            feature[i] = CTDT(sequence)
        if var[i] == 10:
            feature[i] = GTPC(sequence)
        if var[i] == 11:
            feature[i] = TAAI(sequence, use_pro = [87, 379, 29, 88, 145, 30, 25, 26, 27, 95, 452, 92, 154, 15, 81])
        if var[i] == 12:
            feature[i] = GAAC(sequence)

    for i in range(n):
        if i == 0:
            fused_feature = feature[i]
        else:
            fused_feature = np.concatenate( (fused_feature, feature[i]), axis=1)
            
    
    return fused_feature







