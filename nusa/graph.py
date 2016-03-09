# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
#~ """
#~ The graph module contains functions for 
#~ plot displacements, stress, strains and 
#~ other outputs from analysis.
#~ """
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.collections as collections

class GSpring(lines.Line2D):
    def __init__(self,init,**kwargs):
        lines.Line2D.__init__(self,[],[],**kwargs)
        self.xdata = np.array([0,1/5.0,2/5.,3/5.,4/5.,1])+init
        self.ydata = np.array([0,0,0.05,-0.05,0,0])
        # Set data
        self.set_xdata(self.xdata)
        self.set_ydata(self.ydata)
        #self.set_linestyle("-") #Solid line
        
class GFixedNode(collections.PatchCollection):
    def __init__(self,xcoord,ycoord,**kwargs):
        collections.PatchCollection.__init__(self,[],**kwargs)
        tol = 0.02
        _xdata1 = np.array([xcoord-tol,xcoord,xcoord+tol])
        _ydata1 = np.array([ycoord-tol,ycoord,ycoord-tol])
        _xdata2 = np.array([xcoord-tol,xcoord,xcoord-tol])
        _ydata2 = np.array([ycoord-tol,ycoord,ycoord+tol])
        # Polygons
        p1 = patches.Polygon(zip(_xdata1,_ydata1))
        p1.set_color("r")
        p2 = patches.Polygon(zip(_xdata2,_ydata2))
        #print p1,p2
        #p2.set_color("g")
        # Set data
        self.set_paths((p1,p2))
        self.set_color("k")
        #self.set_marker("-")
        #self.set_mfc('r')
        #self.set_ms(10)
        
        
def plot_spring_model(model):
    GRAPH_TOL = 0.1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for k in model.elements.keys():
        ax.add_line(GSpring(k))
    
    #~ # BCs
    for k,v in model.U.items():
        if v==0:
            ax.add_collection(GFixedNode(k,0))
    
    max_x = len(model.elements)
    max_y = 0.2
    ax.set_xlim(-GRAPH_TOL, max_x+GRAPH_TOL)
    ax.set_ylim(-max_y,max_y)
    fig.savefig('this.png')
    plt.show()
    
if __name__ == '__main__':
    pass
