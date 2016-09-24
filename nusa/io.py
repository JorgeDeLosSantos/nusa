# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.blogspot.mx
#  License: MIT License
# ***********************************
"""
The purpose of this module is to provide tools to build 
a model automatically from text files with coordinates 
and connectivities.
"""
import numpy as np
from core import *
from element import *
from model import *


def ModelFromFiles(nodesfile,elementsfile,model):
    """
    Creates a model from ASCII files, where nodesfile contains 
    the coordinates X/Y and elementsfile contains the connectivity 
    of elements.
    """
    pass
    #~ dlm = "," # CSV values
    #~ NC = np.loadtxt(nodesfile, delimiter=dlm)
    #~ EC = np.loadtxt(elementsfile, delimiter=dlm)
    
    #~ for nd in NC:
        #~ cnode = Node( (nd[0], nd[1]) )
        #~ model.addNode(cnode)
        
    #~ for el in EC:
        #~ na = model.nodes[el[0]] 
        #~ nb = model.nodes[el[1]] 
        #~ cel = Spring((na,nb))
        #~ model.addElement(cel)
    
    
if __name__=='__main__':
    pass
