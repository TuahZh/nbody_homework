#########################################################################
# File Name: plot.py
# Author: ZHANG, GaoYuan
# mail: zgy0106@gmail.com
# Created Time: Tue 17 Nov 2015 10:33:36 PM CST
#########################################################################
#description

#!/bin/env python

from __future__ import division
from pylab import *

tmp = loadtxt("2.dat")
data1 = tmp[tmp[:,0]==1]
data2 = tmp[tmp[:,0]==2]
data3 = tmp[tmp[:,0]==3]

print data1.shape, data2.shape, data3.shape

fig = figure(14)
ax = fig.add_subplot(111)
#ax.plot(data1[0:100,1],data1[0:100,2],'r*',data2[0:100,1],data2[0:100,2],'g-', data3[0:100,1],data3[0:100,2],'b-')
ax.plot(data1[0,1],data1[0,2],'r*',data2[0,1],data2[0,2],'g*', data3[0,1],data3[0,2],'b*')
ax.plot(data1[:,1],data1[:,2],'r-',data2[:,1],data2[:,2],'g-', data3[:,1],data3[:,2],'b-')
ax.set_aspect("equal")
#ax.plot(data2[:,1],data2[:,2],'g', data3[:,1],data3[:,2],'b')
fig.savefig("2.png")
