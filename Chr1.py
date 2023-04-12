import csv
import pandas as pd
import screed 
import statistics
import math as m
import numpy as np
import random
import matplotlib.pyplot as plt
import Background
from typing import List, Set
from collections import defaultdict
from collections import OrderedDict

#Dont include anything in 1st row of excel

seqs = Background.readFASTA("C:/Users/Chintu Tinku/Desktop/Project/ChrI SC.fasta")
path= "C:/Users/Chintu Tinku/Desktop/Project/Genes  RefSeq propagation from SGD  annotation version R64-3-1.CSV"
d= pd.read_csv(path)

#Using Hotspot and Coldspot data from Gerton et.al.,2000
path_dataset= "C:/Users/Chintu Tinku/Desktop/Project/Dataset_final.fasta"
recsites = Background.findhotcold(path_dataset)
hot,cold = Background.hotandcold(recsites)
newhotspots,newcoldspots= Background.trimlength(hot,cold)
# ReadFASTA -- detects multiple FASTA files and returns as a single string.
# recsites -- returns a list of untrimmed hot and coldspots sequences
# Newhot/coldspots -- has only the sequence as string


#tetramer analysis

k=4
st=[]
st= Background.allcombinations(st) # to generate all 256 combinations of kmers
kmerhash= Background.kmerhashvalues(st)
hotspots_kmerfreq= Background.matrixforhotspots(newhotspots)
coldspots_kmerfreq= Background.matrixforcoldspots(newcoldspots)

# hotspots_kmerfreq -- 478 x 256 2D array having kmer frequencies of hotspot sequences.
# coldspots_kmerfreq -- 572 x 256 2D array having kmer frequencies of coldspot sequences.


newhotspotlengths=[]
newcoldspotslengths=[]
#A list of lengths of trimmed hotspots and coldspots
for nh in newhotspots:
    newhotspotlengths.append(len(nh))
for nc in newcoldspots:
    newcoldspotslengths.append(len(nc))

trainhotspots=[]
testhotspots=[]
traincoldspots=[]
testcoldspots=[]

hotspotindset = [[429, 453, 1, 254, 16, 19, 422, 276, 73, 205, 20, 447, 435, 425, 374, 267, 82, 419, 70, 50, 297, 304, 284, 306, 237, 132, 313, 206, 400, 252, 74, 178, 474, 145, 312, 83, 131, 7, 94, 454, 28, 41, 369, 244, 264, 48, 194, 340, 322, 18, 461, 335, 176, 346, 434, 442, 49, 242, 214, 290, 427, 197, 345, 165, 381, 98, 426, 281, 424, 438, 470, 90, 328, 289, 212, 384, 268, 401, 292, 311, 152, 26, 446, 192, 116, 353, 34, 5, 220, 53, 162, 137, 76, 36, 21, 227], [158, 410, 150, 368, 2, 354, 279, 231, 318, 79, 89, 258, 370, 31, 25, 259, 161, 173, 416, 323, 316, 121, 160, 451, 199, 393, 330, 130, 472, 469, 184, 198, 201, 189, 402, 255, 275, 344, 203, 249, 343, 66, 159, 56, 24, 385, 392, 241, 319, 163, 271, 367, 151, 371, 358, 327, 432, 236, 187, 123, 12, 113, 200, 364, 127, 174, 40, 72, 195, 133, 77, 372, 277, 101, 307, 310, 78, 437, 280, 406, 321, 439, 303, 148, 54, 190, 
128, 412, 147, 226, 140, 27, 363, 110, 141, 475], [294, 37, 298, 455, 460, 81, 398, 10, 229, 334, 61, 468, 449, 305, 224, 296, 47, 314, 138, 13, 
112, 186, 167, 233, 93, 114, 170, 55, 462, 181, 8, 104, 6, 129, 134, 119, 263, 142, 403, 465, 444, 407, 440, 246, 45, 430, 389, 337, 420, 188, 342, 209, 458, 193, 408, 373, 286, 378, 52, 172, 245, 211, 59, 269, 394, 397, 168, 459, 39, 336, 399, 64, 99, 456, 436, 266, 63, 390, 315, 387, 124, 457, 240, 185, 386, 102, 100, 388, 154, 155, 182, 332, 179, 71, 67, 376], [441, 356, 144, 204, 235, 450, 391, 42, 234, 96, 251, 477, 375, 136, 
51, 473, 308, 239, 295, 404, 417, 84, 32, 4, 414, 396, 15, 278, 62, 293, 60, 208, 395, 366, 146, 285, 3, 180, 380, 452, 232, 105, 283, 409, 260, 
443, 405, 95, 118, 467, 125, 411, 157, 355, 270, 431, 274, 287, 58, 103, 256, 91, 43, 207, 175, 448, 183, 9, 143, 17, 166, 382, 44, 169, 117, 415, 46, 191, 464, 272, 248, 57, 418, 329, 288, 291, 111, 115, 218, 202, 309, 177, 261, 250, 349, 35], [0, 11, 14, 22, 23, 29, 30, 33, 38, 65, 68, 69, 75, 80, 85, 86, 87, 88, 92, 97, 106, 107, 108, 109, 120, 122, 126, 135, 139, 149, 153, 156, 164, 171, 196, 210, 213, 215, 216, 217, 219, 221, 
222, 223, 225, 228, 230, 238, 243, 247, 253, 257, 262, 265, 273, 282, 299, 300, 301, 302, 317, 320, 324, 325, 326, 331, 333, 338, 339, 341, 347, 
348, 350, 351, 352, 357, 359, 360, 361, 362, 365, 377, 379, 383, 413, 421, 423, 428, 433, 445, 463, 466, 471, 476]]

coldspotindset = [[156, 418, 154, 477, 39, 85, 302, 204, 165, 35, 513, 504, 485, 360, 421, 226, 420, 380, 86, 110, 111, 313, 269, 437, 411, 299, 456, 257, 113, 497, 470, 423, 222, 484, 163, 73, 489, 363, 529, 438, 69, 115, 239, 544, 567, 57, 550, 264, 406, 322, 275, 80, 417, 122, 11, 384, 127, 287, 55, 143, 471, 147, 242, 245, 546, 352, 286, 448, 174, 46, 327, 337, 453, 232, 526, 219, 510, 179, 48, 197, 157, 543, 432, 458, 117, 534, 176, 549, 316, 
60, 542, 188, 95, 230, 276, 416, 521, 206, 141, 203, 148, 317, 58, 558, 30, 37, 433, 425, 59, 142, 32, 412, 463, 17], [507, 212, 271, 435, 440, 210, 295, 205, 255, 216, 145, 350, 248, 62, 121, 353, 14, 340, 246, 220, 441, 218, 92, 131, 311, 288, 318, 405, 81, 90, 189, 134, 486, 207, 398, 293, 27, 495, 379, 152, 266, 312, 136, 166, 422, 247, 300, 159, 71, 240, 227, 158, 375, 31, 231, 249, 201, 532, 474, 7, 364, 114, 233, 407, 224, 386, 265, 480, 445, 490, 344, 496, 193, 194, 252, 568, 65, 20, 525, 202, 482, 263, 182, 97, 368, 25, 511, 554, 191, 137, 324, 292, 254, 91, 63, 54, 332, 236, 476, 297, 12, 23, 184, 478, 301, 267, 244, 116, 50, 395, 414, 164, 126, 229], [74, 49, 211, 498, 381, 162, 540, 234, 535, 34, 33, 303, 284, 494, 268, 462, 372, 6, 161, 29, 348, 570, 444, 135, 520, 394, 431, 424, 225, 281, 15, 120, 124, 564, 178, 449, 170, 105, 79, 413, 370, 566, 354, 153, 18, 427, 541, 3, 319, 102, 8, 377, 93, 409, 66, 436, 183, 461, 467, 509, 404, 367, 2, 132, 198, 171, 347, 447, 376, 101, 44, 175, 290, 190, 51, 338, 468, 291, 483, 410, 72, 64, 294, 192, 460, 77, 465, 351, 9, 96, 388, 342, 357, 419, 397, 374, 539, 106, 362, 524, 305, 237, 538, 
334, 531, 108, 213, 488, 4, 528, 243, 167, 168, 235], [343, 396, 553, 333, 314, 56, 70, 238, 563, 450, 41, 457, 464, 402, 38, 446, 329, 451, 130, 253, 87, 149, 371, 527, 36, 277, 551, 285, 160, 282, 274, 82, 328, 100, 366, 392, 181, 112, 98, 272, 335, 555, 298, 5, 506, 200, 308, 0, 52, 561, 501, 283, 492, 349, 256, 261, 22, 45, 400, 1, 365, 251, 10, 83, 514, 104, 84, 273, 475, 339, 361, 40, 569, 505, 13, 455, 138, 250, 61, 479, 500, 330, 442, 517, 150, 139, 454, 199, 129, 434, 383, 270, 481, 78, 321, 512, 47, 515, 439, 128, 571, 144, 491, 552, 228, 523, 208, 241, 99, 499, 186, 562, 356, 530], [16, 19, 21, 24, 26, 28, 42, 43, 53, 67, 68, 75, 76, 88, 89, 94, 103, 107, 109, 118, 119, 123, 125, 133, 140, 146, 151, 155, 
169, 172, 173, 177, 180, 185, 187, 195, 196, 209, 214, 215, 217, 221, 223, 258, 259, 260, 262, 278, 279, 280, 289, 296, 304, 306, 307, 309, 310, 
315, 320, 323, 325, 326, 331, 336, 341, 345, 346, 355, 358, 359, 369, 373, 378, 382, 385, 387, 389, 390, 391, 393, 399, 401, 403, 408, 415, 426, 
428, 429, 430, 443, 452, 459, 466, 469, 472, 473, 487, 493, 502, 503, 508, 516, 518, 519, 522, 533, 536, 537, 545, 547, 548, 556, 557, 559, 560, 
565]]


hotspotIDset = []
coldspotIDset = []
lhset = []
indtesthset = []
lcset = []
indtestcset = []


for valloop in range(0,5):
    lh = []     #lh-> indexes of training data hotspots-> 382 hotspots
    lc = []     #lc-> indexes of training data coldspots-> 458 coldspots
    indtesth = []
    indtestc = []
    for i in range(0,5):
        if i == valloop: 
            indtesth+=hotspotindset[i]
            indtestc+=coldspotindset[i]
            continue
        else:
            lh+= hotspotindset[i]
            lc+= coldspotindset[i]
    lhset.append(lh)
    lcset.append(lc)
    indtesthset.append(indtesth)
    indtestcset.append(indtestc)

    #----------------------------Increment of Diversity----------------------------#
    reptrainh= [0]*256
    reptrainc= [0]*256

    for i in range(0,256):
        for j in lh:
            reptrainh[i]+= hotspots_kmerfreq[j][i]

    for i in range(0,256):
        for j in lc:
            reptrainc[i]+= coldspots_kmerfreq[j][i]

    #test = coldspots_kmerfreq[470]

    allID=[] #List of 2D tuple with ID value between hotspot training and coldspot training class 
    hotID=[]
    coldID=[]
    kmerfreqtotal = hotspots_kmerfreq + coldspots_kmerfreq

    for i in kmerfreqtotal:
        allID.append((Background.IncDiv(i,reptrainh),Background.IncDiv(i,reptrainc)))

    for i in hotspots_kmerfreq:
        hotID.append((Background.IncDiv(i,reptrainh),Background.IncDiv(i,reptrainc)))

    for i in coldspots_kmerfreq:
        coldID.append((Background.IncDiv(i,reptrainh),Background.IncDiv(i,reptrainc)))

    hotspotIDset.append(hotID)
    coldspotIDset.append(coldID)


'''
#---------------------Quadratic Discrimination Analysis-------------------------#

#lh ---index of training class hotspots
#lc ---index of training class coldspots

hotIDtrain =[]
coldIDtrain =[]

for i in lh:
    hotIDtrain.append(hotID[i])
for i in lc:
    coldIDtrain.append(coldID[i])

# print(len(hotIDtrain)) --> 382 sequence
# print(len(coldIDtrain)) --> 458 sequence

testID = coldID[17] # 470 is the index of a test sequence in coldspot class

meantrainIDhot = Background.meanIDvctr(hotIDtrain)
#print(indtesth,indtestc)
meantrainIDcold = Background.meanIDvctr(coldIDtrain)

##------Covariance matrix generation--------
  # Generate the matrix seperately for hotspot and coldspot classes seperately

covIDhp = [[0]*2 for _ in range(2)]
covIDcp = [[0]*2 for _ in range(2)]

for j in range(0,2):
    for i in range(0,2):
        value=0
        for s in lh:
            value+= ((hotID[s][j]- meantrainIDhot[j]) *(hotID[s][i]- meantrainIDhot[i]))
        covIDhp[j][i]= (value/(len(lh)-1))

for j in range(0,2):
    for i in range(0,2):
        value=0
        for s in lc:
            value+= ((coldID[s][j]- meantrainIDcold[j]) *(coldID[s][i]- meantrainIDcold[i]))
        covIDcp[j][i]= (value/(len(lc)-1))

##------Mehalanobis distance--------

mehdisth = Background.mehdist(testID, meantrainIDhot, covIDhp)
mehdistc = Background.mehdist(testID, meantrainIDcold, covIDcp)

##------Finding the Quadratic Discrimination function------

lenhtrain = len(hotIDtrain)
lenctrain = len(coldIDtrain)

detcovIDhp = Background.determinant(covIDhp)
detcovIDcp = Background.determinant(covIDcp)

QDF = (m.log(lenhtrain/lenctrain,2))-((mehdisth-mehdistc)/2)-((m.log(detcovIDhp/detcovIDcp,2))/2)

#print(QDF)

kmerfreqcoldall= [0]*256 # kmer frequenzy of all coldspot training class sequences.
for i in range(0,256):
    value=0
    for s in coldspots_kmerfreq:
        value+= s[i]
    kmerfreqcoldall[i]= value
#print(kmerfreqcoldall)'''



