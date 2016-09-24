#~ # -*- coding: utf-8 -*-
#~ # ***********************************
#~ #  Author: Pedro Jorge De Los Santos     
#~ #  E-mail: delossantosmfq@gmail.com 
#~ #  Web: labdls.blogspot.mx
#~ #  License: MIT License
#~ # ***********************************
#~ import sys
#~ sys.path.insert(0, '../')
#~ 
#~ import matplotlib.pyplot as plt
from nusa.core import *
#~ from nusa.lib import STEEL_1018, ALUMINIUM_6061
from nusa.model import SpringModel
from nusa.element import Spring
#~ from nusa.graph import plot_spring_model

import nusa.io as io

m1 = SpringModel()
a = io.ModelFromFiles("data/nodes.txt","data/elements.txt", m1)

for nd in m1.getNodes():
    print nd.x, nd.y
