#! /usr/bin/env python
#!_*_ coding: utf-8

#==========================
# @2017.01.12
# @hzau.wuhan.china
# @liufang
#==========================


import os
import numpy as np


from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs


data_file = "a.out"
X = []
index = []
fp = open(data_file, 'r')
line = fp.readline().strip()
line = fp.readline().strip()
while line != "":
    line = line.split()
    index.append(line[0])
    rec = [float(r) for r in line[2:] ]
    X.append(rec)
    line = fp.readline().strip()
fp.close()
X = np.array(X)

y_pred = KMeans(n_clusters=12).fit_predict(X)

out_file = 'a.dat'
fp = open(out_file, 'w')
for i in range(len(index)):
    x = ["%8s" %str(m) for m in X[i]]
    x = "  ".join(x)
    print >>fp, "%-10s%s%10d%1s" %(index[i],x, y_pred[i],'')
fp.close()


