# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import numpy as np

nodos = np.loadtxt("nlist")[:,1:3]
elementos = np.loadtxt("elist")[:,6:9]
np.savetxt("nodos",nodos)
np.savetxt("elementos",elementos,fmt="%i")
