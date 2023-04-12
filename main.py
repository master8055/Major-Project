import csv
import Background
import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics as stat
import chr1


AccListset = []
xROCset = []
yROCset = []

for valloop in range(0,5):

    hotID = chr1.hotspotIDset[valloop]
    coldID = chr1.coldspotIDset[valloop]

    lh = chr1.lhset[valloop]
    lc = chr1.lcset[valloop]

    indtesth = chr1.indtesthset[valloop]
    indtestc = chr1.indtestcset[valloop]



    #---------------------Quadratic Discrimination Analysis-------------------------#

    #lh ---index of training class hotspots
    #lc ---index of training class coldspots

    hotIDtrain =[]
    coldIDtrain =[]

    for i in lh:
        hotIDtrain.append(hotID[i])
    for i in lc:
        coldIDtrain.append(coldID[i])

    meantrainIDhot = Background.meanIDvctr(hotIDtrain)
    meantrainIDcold = Background.meanIDvctr(coldIDtrain)

    # print(len(hotIDtrain)) --> 382 sequence
    # print(len(coldIDtrain)) --> 458 sequence

    # List of 2D tuple with (index, QDF) for test sequences
    QDFH = []
    QDFC = []

    QDFhot=[] # Has only the QDF value for hotspot test sequences
    QDFcold=[] # Has only the QDF value for coldspot test sequences

    #Mahalanobis distance of hotspot test sequences with hotspot and coldspot training class respectively
    mdistlisthh=[]
    mdistlisthc=[]
    #Mahalanobis distance of coldspot test sequences with hotspot and coldspot training class respectively
    mdistlistch=[]
    mdistlistcc=[]


    ##------Covariance matrix generation--------
    #-------Generate the matrix seperately for hotspot and coldspot training classes seperately

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

    #--------------------------------------------

    for k in indtesth:

        testID = hotID[k]
        
        ##------Mehalanobis distance--------

        mehdisth = Background.mehdist(testID, meantrainIDhot, covIDhp)
        mehdistc = Background.mehdist(testID, meantrainIDcold, covIDcp)

        mdistlisthh.append(mehdisth)
        mdistlisthc.append(mehdistc)
        ##------Finding the Quadratic Discrimination function------

        lenhtrain = len(hotIDtrain)
        lenctrain = len(coldIDtrain)

        detcovIDhp = Background.determinant(covIDhp)
        detcovIDcp = Background.determinant(covIDcp)

        QDF = (m.log(lenhtrain/lenctrain,2))-((mehdisth-mehdistc)/2)-((m.log(detcovIDhp/detcovIDcp,2))/2)
        QDFhot.append(QDF)
        QDFH.append((k,QDF))

    for k in indtestc:

        testID = coldID[k]

        ##------Mehalanobis distance--------

        mehdisth = Background.mehdist(testID, meantrainIDhot, covIDhp)
        mehdistc = Background.mehdist(testID, meantrainIDcold, covIDcp)
        mdistlistch.append(mehdisth)
        mdistlistcc.append(mehdistc)

        ##------Finding the Quadratic Discrimination function------

        lenhtrain = len(hotIDtrain)
        lenctrain = len(coldIDtrain)

        detcovIDhp = Background.determinant(covIDhp)
        detcovIDcp = Background.determinant(covIDcp)

        QDF = (m.log(lenhtrain/lenctrain,2))-((mehdisth-mehdistc)/2)-((m.log(detcovIDhp/detcovIDcp,2))/2)
        QDFcold.append(QDF)
        QDFC.append((k,QDF))
    #QDF value for all hotspot and coldspot test sequences are generated and stored in QDFH and QDFC respectively.

    #-----------------------------Statistics-----------------------------------#

    xROC =[]
    yROC=[]
    AccList=[]
    Ethrlist = [i/10 for i in range(-100,40,1)]
    for Ethr in Ethrlist:
        TP = 0
        TN = 0
        FP = 0
        FN = 0
        for i in QDFH:
            if i[1]>Ethr:
                TP+=1
            else:
                FN+=1
        for i in QDFC:
            if i[1]<=Ethr:
                TN+=1
            else:
                FP+=1

        Sn= (TP)/(TP+FN)
        Sp= (TN)/(TN+FP)
        TA= (TP+TN)/(TP+TN+FP+FN)
        a= (TP*TN)-(FP*FN)
        b= m.sqrt((TP+FP)*(TN+FN)*(TP+FN)*(TN+FP))
        #CC= a/b
        xROC.append(1-Sp)
        yROC.append(Sn)
        AccList.append(TA*100)

    xROCset.append(xROC)
    yROCset.append(yROC)
    AccListset.append(AccList)





'''for i in range(0,5):
    plt.plot(xROCset[i], yROCset[i], label = 'set_'+ str(i+1))

plt.xlabel("1-Specificity(1-Sp)")
plt.ylabel("Sensitivity(Sn)")
plt.title("Receiver Operating Characteristics(ROC) of all 5 sets")
plt.legend()
plt.show()
'''

'''
#plt.yticks([stat.median(QDFH),5,10,15,20,25,30])
plt.ylabel('Quadratic Discrimination Function(QDF)')  
plt.title('The boxplot of predicted QDF value of 114 coldspot test sequences ')
plt.boxplot(QDFcold)
plt.yticks([max(QDFcold),min(QDFcold),stat.median(QDFcold),-12,-10,-8,-6,-4,-2,0,2,4,6])
plt.show()'''




