import main
import pandas as pd
import screed 
import numpy as np
import math as m
import pandas as pd
import random

hotspotset = []
coldspotset = []
hotspotset.append(main.indtesth)
coldspotset.append(main.indtestc)

lh= main.lh
indtesth= []
testhotspots=[]

for j in range(3):
    indtesth= []
    for i in range(1,97):
        index= random.choice(lh)
        lh.remove(index)
        indtesth.append(index)
            
    hotspotset.append(indtesth)

hotspotset.append(lh)
lc= main.lc

for j in range(3):
    indtestc= []
    for i in range(1,115):
        index= random.choice(lc)
        lc.remove(index)
        indtestc.append(index)
            
    coldspotset.append(indtestc)

coldspotset.append(lc)

for i in coldspotset:
    print(len(i))





