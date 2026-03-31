#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 11:11:44 2021

@author: tatianamiti
"""


####### gives a random distribution of cells given a tissue sample with stroma, positivie and negative brdu cells. Chooses randomly brdu+ cells from all cells,
########### the rest are brdu-

import os,json
import csv
import random, math
import matplotlib.pyplot as plt
import shapely

if os.path.isfile('numberFile.txt'):
    with open('numberFile.txt', 'r') as f_f:
        numberFile = json.loads(f_f.read())
else:
    numberFile = []

numb = numberFile[0]
print(numb)
    
if os.path.isfile('StromaList.txt'):
    with open('StromaList.txt', 'r') as fff:
        StromaList = json.loads(fff.read())
else:
    StromaList = []
    

if os.path.isfile('NegativeList.txt'):
        with open('NegativeList.txt', 'r') as f_ff:
            NegativeList = json.loads(f_ff.read())
else:
            NegativeList = []
            
            
if os.path.isfile('PositiveList.txt'):
        with open('PositiveList.txt', 'r') as f_f:
            PositiveList = json.loads(f_f.read())
else:
            PositiveList = []
    
print((len(PositiveList) + len(NegativeList) + len(StromaList)))
    


## create lists to keep the data in separate lists that will be individuallu modified
ListAll = []
print("pos and neg inside a tissue piece", len(PositiveList), len(NegativeList))

### list used to assign labels
for r in PositiveList:
    ListAll.append([r[0], r[1]])
    
for y in NegativeList:
    ListAll.append([y[0], y[1]])
### label brdu+ cells within 2 diametersfrom stroma, brdu_ the rest
tempPos = []
tempNeg = []
    
    
while len(tempPos) <= len(PositiveList):
    b = random.choice(ListAll)
    if b not in tempPos:
        tempPos.append(b)
        ListAll.remove(b)
        
print("got the positives", len(tempPos))
tempNeg = ListAll[:]
# for rr in ListAll:
#     if rr not in tempPos:
#         tempNeg.append(rr)
    
print("positive new, old", len(tempPos), len(PositiveList), len(tempNeg), len(NegativeList))
ListAll = []
Q1 = []
for nm in tempNeg:
    Q1.append(['BrdU_neg',nm[0], nm[1]])
for mk in StromaList:
    Q1.append(['Stroma', mk[0], mk[1]])
for bn in tempPos:
    Q1.append(['BrdU_pos', bn[0], bn[1]])
          
    
print("pos neg str all cells =", len(tempPos),len(tempNeg), len(StromaList), len(Q1))
    # write the listst to csv files for R
labels = ['class_label','X', 'Y']

with open(str(numb) + '_Converted_coordinates_DIAMETER_TRUERand.csv', 'w') as f: 
    write = csv.writer(f) 
    write.writerow(labels) 
    write.writerows(Q1) 
    