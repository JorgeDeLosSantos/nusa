# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.tri as tri

fig = plt.figure()
ax = fig.add_subplot(111)

x,y = [0,1,0],[0,0.5,1.0]
vtx = (20,20,10)

plt.subplot(111)
triang = tri.Triangulation(x, y)
print dir(triang)
plt.tricontourf(x, y, vtx, cmap="jet")

plt.show()
