# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 21:24:48 2016
use python 3
"""

import csv

compare = [3,5,7]

f = open('resultcsv.csv', 'r')
reader = csv.reader(f)
header = next(reader)
tmp = 0
indx = 0
min_err = 1.0E-1
ans= {1:0,2:0,3:0,4:0,5:0,6:0}
count = 0
for row in reader:
    if(row[0].find("max")!=-1):
        print(row[0])
    else:
        method = int(row[1])
        err = float(row[3])
        last = min(compare)
        if method in compare :
            if(min_err > err):
                indx = method
                min_err = err
            if(method == last):
                print(indx,min_err)
                ans[indx] += 1
                count += 1
                indx = 0
                min_err = 10.0


print("割合は...(",count,"個のインスタンス中)")
for x in ans:
    if ans[x]!=0:
        print(x,ans[x]/count)
f.close()
