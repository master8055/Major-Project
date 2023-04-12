from typing import List, Set
import csv
import pandas as pd
import screed 
import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import OrderedDict

recsites=[]
hotspots=[]
coldspots=[]
sun: int

#------------------------Processing---------------------
def readFASTA(inputfile): 
    with open(inputfile, "r") as f: 
        seq = f.readline()
        seq = f.read() 
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
        
    return seq 

def findhotcold(path):

    seq= readFASTA(path)
    cnt=0
    recsites=[]
    st=""
    for i in seq:
        if(cnt == 0):
            cnt+=1
            continue

        if(i != '>'):
            st+= i

        if(i == '>'):
            recsites.append(st)
            st=""

    recsites.append(st)
    return recsites

def hotandcold(recsites):
    for i in range(0,len(recsites)):
        if(recsites[i].find("Hot") != -1):
            hotspots.append(recsites[i])
        else:
            coldspots.append(recsites[i])   

    return hotspots,coldspots

def trimlength(hot, cold):
    index=0
    newhotspots=[]
    newcoldspots=[]
    for i in range(0, len(hot)):
        index= hot[i].find("Hot")
        index+= 3
        x= hot[i]
        newhotspots.append(x[index:])

    for i in range(0, len(cold)):
        index= cold[i].find("Cold")
        index+= 4
        x= cold[i]
        newcoldspots.append(x[index:])
    return newhotspots,newcoldspots

#---------------------------------------------------------

s=[]
v: List[str] = ['A', 'C', 'G', 'T']
kmerhash = OrderedDict()


#--Generating all possible k-mers and hashing it from 0 to 255--

def generate(i: int, ans: str, s) -> None:
    if i == 4:
        s.append(ans)
        return
    for it in v:
        generate(i+1, ans+it, s)

def allcombinations(st):
    generate(0, "",s)
    s.sort()
    st= s
    return st

def kmerhashvalues(s):
    start = 0
    for it in s:
        kmerhash[it] = start
        start += 1
    #for it in kmerhash.items():
        #print(it[0], it[1])
    return kmerhash

#-----------------------------------------------------------------


#--Generate k-mer frequencies of all hotspot and coldspot sequencies--

def matrixforhotspots(newhotspots):
    ind = 0
    k=4
    hotspots_kmerfreq = [[0]*256 for _ in range(478)]
    for i in newhotspots:
        for j in range(0,len(i)-k+1):
            seq= i[j:j+k] #generated kmer for ith hotspot
            hotspots_kmerfreq[ind][kmerhash[seq]]+= 1
        
        ind+= 1
    return hotspots_kmerfreq

def matrixforcoldspots(newcoldspots):
    ind=0
    k=4
    coldspots_kmerfreq = [[0]*256 for _ in range(572)]
    for i in newcoldspots:
        for j in range(0,len(i)-k+1):
            seq= i[j:j+k]
            coldspots_kmerfreq[ind][kmerhash[seq]]+= 1
            #generated kmer for ith hotspot
        ind+= 1
    return coldspots_kmerfreq

#------------------------------------------------------------------



def multiply(a, b):
    return np.dot(a,b)

def add(a, b):
    return np.add(a,b)

def transpose(a):
    return np.transpose(a)

def inverse(a):
    return np.linalg.inv(a)

def determinant(a):
    return np.linalg.det(a)

def subtract(a,b):
    return np.subtract(a,b)


def diversity(sequence:List):
    divN= 0
    for i in sequence:
        divN+= i
    divA=0
    divB=0
    divA= divN*(np.log(divN))

    for i in sequence:
        if i !=0:
            divB+= i*(np.log(i))
        else: 
            continue
    return divA-divB

def IncDiv(A:List, B:List):
    C=[0]*256
    for i in range(0,255):
        C[i] = A[i]+B[i]

    IncDiv= diversity(C) - diversity(A) - diversity(B)
    return IncDiv

def meanIDvctr(X: List):

    sumh=0
    sumc=0
    meanh = 0
    meanc = 0
    length = len(X)
    for i in X:
        sumh += i[0]
        sumc += i[1]
    meanh = sumh/length
    meanc = sumc/length

    return [meanh,meanc]

def mehdist(test, meanIDvctr, covmat):

    A = transpose(subtract(test,meanIDvctr))
    B = inverse(covmat)
    D = subtract(test,meanIDvctr)

    C = multiply(A,B)

    return multiply(C,D)
    



"""covhp= [[0]*256 for _ in range(256)] #covariance matrix for hotspots
covcp= [[0]*256 for _ in range(256)] #covariance matrix for coldspots

for j in range(0,256):
    for i in range(0,256):
        value=0
        for s in lh:
            value+= ((hotspots_kmerfreq[s][j]- meanvectorhp[j]) *(hotspots_kmerfreq[s][i]- meanvectorhp[i]))
        covhp[j][i]= (value/(len(lh)-1))

for j in range(0,256):
    for i in range(0,256):
        value=0
        for s in lc:
            value+= ((coldspots_kmerfreq[s][j]- meanvectorcp[j]) *(coldspots_kmerfreq[s][i]- meanvectorcp[i]))
        covcp[j][i]= value/(len(lc)-1)

#print(np.amax(covhp))
"""

'''l= [0, 8, 128, 32, 3, 2, 64, 63, 255]
kmers=[]
for j in l:
    value = {i for i in kmerhash if kmerhash[i] == j} # kmerhash - ordered dictionary of all kmers
    print(value)
'''

'''
lh= [x for x in range(0,478)] #lh-> indexes of training data hotspots
for i in range(1,97):
    index= random.choice(lh)
    if(index not in indtesth):
        testhotspots.append(newhotspots[index])
        indtesth.append(index)
        lh.remove(index)

lc= [x for x in range(0,572)] #lc-> indexes of training data coldspots
for i in range(1,115):
    index= random.choice(lc)
    if(index not in indtestc):
        testcoldspots.append(newcoldspots[index])
        indtestc.append(index)
        lc.remove(index)
'''


            

