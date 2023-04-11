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


#------------tetramer analysis-------------#

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
indtesth= [429, 453, 1, 254, 16, 19, 422, 276, 73, 205, 20, 447, 435, 425, 374, 267, 82, 419, 70, 50, 297, 304, 284, 306, 237, 132, 313, 206, 400, 252, 74, 178, 474, 145, 312, 83, 131, 7, 94, 454, 28, 41, 369, 244, 264, 48, 194, 340, 322, 18, 461, 335, 176, 346, 434, 442, 49, 242, 214, 290, 427, 197, 345, 165, 381, 98, 426, 281, 424, 438, 470, 90, 328, 289, 212, 384, 268, 401, 292, 311, 152, 26, 446, 192, 116, 353, 34, 5, 220, 53, 162, 137, 76, 36, 21, 227]
# index of test hotspots 96
indtestc= [156, 418, 154, 477, 39, 85, 302, 204, 165, 35, 513, 504, 485, 360, 421, 226, 420, 380, 86, 110, 111, 313, 269, 437, 411, 299, 456, 257, 113, 497, 470, 423, 222, 484, 163, 73, 489, 363, 529, 438, 69, 115, 239, 544, 567, 57, 550, 264, 406, 322, 275, 80, 417, 122, 11, 384, 127, 287, 55, 143, 471, 147, 242, 245, 546, 352, 286, 448, 174, 46, 327, 337, 453, 232, 526, 219, 510, 179, 48, 197, 157, 543, 432, 458, 117, 534, 176, 549, 316, 60, 542, 188, 95, 230, 276, 416, 521, 206, 141, 203, 148, 317, 58, 558, 30, 37, 433, 425, 59, 142, 32, 412, 463, 17]
# index of test coldspots 114

lh = [0, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 17, 22, 23, 24, 25, 27, 29, 30, 31, 32, 33, 35, 37, 38, 39, 40, 42, 43, 44, 45, 46, 47, 51, 52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 71, 72, 75, 77, 78, 79, 80, 81, 84, 85, 86, 87, 88, 89, 91, 92, 
93, 95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 117, 118, 119, 120, 121, 122, 123, 
124, 125, 126, 127, 128, 129, 130, 133, 134, 135, 136, 138, 139, 140, 141, 142, 143, 144, 146, 147, 148, 149, 150, 151, 153, 154, 155, 
156, 157, 158, 159, 160, 161, 163, 164, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 177, 179, 180, 181, 182, 183, 184, 185, 186, 
187, 188, 189, 190, 191, 193, 195, 196, 198, 199, 200, 201, 202, 203, 204, 207, 208, 209, 210, 211, 213, 215, 216, 217, 218, 219, 221, 
222, 223, 224, 225, 226, 228, 229, 230, 231, 232, 233, 234, 235, 236, 238, 239, 240, 241, 243, 245, 246, 247, 248, 249, 250, 251, 253, 
255, 256, 257, 258, 259, 260, 261, 262, 263, 265, 266, 269, 270, 271, 272, 273, 274, 275, 277, 278, 279, 280, 282, 283, 285, 286, 287, 
288, 291, 293, 294, 295, 296, 298, 299, 300, 301, 302, 303, 305, 307, 308, 309, 310, 314, 315, 316, 317, 318, 319, 320, 321, 323, 324, 
325, 326, 327, 329, 330, 331, 332, 333, 334, 336, 337, 338, 339, 341, 342, 343, 344, 347, 348, 349, 350, 351, 352, 354, 355, 356, 357, 
358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 370, 371, 372, 373, 375, 376, 377, 378, 379, 380, 382, 383, 385, 386, 387, 388, 
389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 
418, 420, 421, 423, 428, 430, 431, 432, 433, 436, 437, 439, 440, 441, 443, 444, 445, 448, 449, 450, 451, 452, 455, 456, 457, 458, 459, 
460, 462, 463, 464, 465, 466, 467, 468, 469, 471, 472, 473, 475, 476, 477]
#lh-> indexes of training data hotspots
lc = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 33, 34, 36, 38, 40, 41, 42, 
43, 44, 45, 47, 49, 50, 51, 52, 53, 54, 56, 61, 62, 63, 64, 65, 66, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78, 79, 81, 82, 83, 84, 87, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 112, 114, 116, 118, 119, 120, 121, 123, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 144, 145, 146, 149, 150, 151, 152, 153, 155, 158, 159, 160, 161, 162, 164, 166, 167, 168, 169, 170, 171, 172, 173, 175, 177, 178, 180, 181, 182, 183, 184, 185, 186, 187, 189, 190, 191, 192, 193, 194, 195, 196, 198, 199, 200, 201, 202, 205, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 220, 221, 223, 224, 225, 227, 228, 229, 231, 233, 234, 235, 236, 237, 238, 240, 241, 243, 244, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 258, 259, 260, 261, 262, 263, 265, 266, 267, 268, 270, 271, 272, 273, 274, 277, 278, 279, 280, 281, 282, 283, 284, 285, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 300, 301, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 314, 315, 318, 319, 320, 321, 323, 324, 325, 326, 328, 329, 330, 331, 332, 333, 334, 335, 336, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 353, 354, 355, 356, 357, 358, 359, 361, 362, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 381, 382, 383, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 407, 408, 409, 410, 413, 414, 415, 419, 422, 424, 426, 427, 428, 429, 430, 431, 434, 435, 436, 439, 440, 441, 442, 443, 444, 445, 446, 447, 449, 450, 451, 452, 454, 455, 457, 459, 460, 461, 462, 464, 465, 466, 467, 468, 469, 472, 473, 474, 475, 476, 478, 479, 480, 481, 482, 483, 486, 487, 488, 490, 491, 492, 493, 494, 495, 496, 498, 499, 500, 501, 502, 503, 505, 506, 507, 508, 509, 511, 512, 514, 515, 516, 517, 518, 519, 520, 522, 523, 524, 525, 527, 528, 530, 531, 532, 533, 535, 536, 537, 538, 539, 540, 541, 545, 547, 548, 551, 552, 553, 554, 555, 556, 557, 559, 560, 561, 562, 563, 564, 565, 566, 568, 569, 570, 571]
#lc-> indexes of training data coldspots



#lh-> index of training data hotspots-> 382 hotspots
#lc-> index of training data coldspots> 458 coldspots


#----------------------------Increment of Diversity----------------------------#
reptrainh= [0]*256
reptrainc= [0]*256

for i in range(0,256):
    for j in lh:
        reptrainh[i]+= hotspots_kmerfreq[j][i]

for i in range(0,256):
    for j in lc:
        reptrainc[i]+= coldspots_kmerfreq[j][i]

test = coldspots_kmerfreq[470]

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
#print(kmerfreqcoldall)



