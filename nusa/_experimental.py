# ***********************************
#  Author: Pedro Jorge De Los Santos    
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import numpy as np
import numpy.linalg as la
import json
from nusa import *

# class NusaModelReader(object):
#   def __init__(self,filename):
#       self.filename = filename

def read_truss_model(filename):
    json_file = filename
    with open(json_file, 'r') as myfile:
        data=myfile.read()
    obj = json.loads(data)

    nodes_data = _get_nodes(obj)
    elements_data = _get_elements(obj)
    constraints_data = _get_constraints(obj)
    forces_data = _get_forces(obj)

    nc = nodes_data
    ec = elements_data
    x,y = nc[:,0], nc[:,1]

    nodes = []
    elements = []

    for k,nd in enumerate(nc):
        cn = Node((x[k],y[k]))
        nodes.append(cn)
        
    for k,elm in enumerate(ec):
        i,j,E,A = int(elm[0]-1),int(elm[1]-1),elm[2],elm[3]
        ni,nj = nodes[i],nodes[j]
        ce = Truss((ni,nj), E, A)
        elements.append(ce)
        
    model = TrussModel("Truss Model")
    for n in nodes: model.add_node(n)
    for e in elements: model.add_element(e)
    
    for c in constraints_data:
        k,ux,uy = int(c[0]),c[1],c[2]
        if ~np.isnan(ux) and ~np.isnan(uy):
            model.add_constraint(nodes[k-1], ux=ux, uy=uy)
        elif ~np.isnan(ux):
            model.add_constraint(nodes[k-1], ux=ux)
        elif ~np.isnan(uy):
            model.add_constraint(nodes[k-1], uy=uy)
    
    for f in forces_data:
        k,fx,fy = int(f[0]),f[1],f[2]
        model.add_force(nodes[k-1],(fx,fy))

    return model



def read_spring_model(filename):
    json_file = filename
    with open(json_file, 'r') as myfile:
        data=myfile.read()
    obj = json.loads(data)
    
    nodes_data = _get_nodes(obj)
    elements_data = _get_elements_spring(obj)
    constraints_data = _get_constraints(obj)
    forces_data = _get_forces(obj)

    nc = nodes_data
    ec = elements_data
    x,y = nc[:,0], nc[:,1]

    nodes = []
    elements = []

    for k,nd in enumerate(nc):
        cn = Node((x[k],y[k]))
        nodes.append(cn)
        
    for k,elm in enumerate(ec):
        i,j,ke = int(elm[0]-1),int(elm[1]-1),elm[2]
        ni,nj = nodes[i],nodes[j]
        ce = Spring((ni,nj), ke)
        elements.append(ce)
        
    model = SpringModel("Truss Model")
    for n in nodes: model.add_node(n)
    for e in elements: model.add_element(e)
    
    for c in constraints_data:
        k,ux,uy = int(c[0]),c[1],c[2]
        if ~np.isnan(ux) and ~np.isnan(uy):
            model.add_constraint(nodes[k-1], ux=ux, uy=uy)
        elif ~np.isnan(ux):
            model.add_constraint(nodes[k-1], ux=ux)
        elif ~np.isnan(uy):
            model.add_constraint(nodes[k-1], uy=uy)
    
    for f in forces_data:
        k,fx,fy = int(f[0]),f[1],f[2]
        model.add_force(nodes[k-1],(fx,fy))

    return model



def _get_nodes(obj):
    nn = len(obj["nodes"])
    nodes = np.zeros((nn,2))
    for i,m in enumerate(obj["nodes"]):
        nodes[i,0] = m["x"]
        nodes[i,1] = m["y"]
    return nodes

def _get_elements(obj):
    nn = len(obj["elements"])
    elements = np.zeros((nn,4))
    for i,m in enumerate(obj["elements"]):
        elements[i,0] = m["ni"]
        elements[i,1] = m["nj"]
        elements[i,2] = m["E"]
        elements[i,3] = m["A"]
    return elements

def _get_elements_spring(obj):
    nn = len(obj["elements"])
    elements = np.zeros((nn,3))
    for i,m in enumerate(obj["elements"]):
        elements[i,0] = m["ni"]
        elements[i,1] = m["nj"]
        elements[i,2] = m["ke"]
    return elements


def _get_constraints(obj):
    nn = len(obj["constraints"])    
    const = np.zeros((nn,3))
    for i,m in enumerate(obj["constraints"]):
        const[i,1] = np.nan if m["ux"] == "free" else m["ux"]
        const[i,2] = np.nan if m["uy"] == "free" else m["uy"]
        const[i,0] = m["node"]
    return const

def _get_forces(obj):
    nn = len(obj["forces"])
    forces = np.zeros((nn,3))
    for i,m in enumerate(obj["forces"]):
        forces[i,1] = m["fx"]
        forces[i,2] = m["fy"]
        forces[i,0] = m["node"]
    return forces


if __name__=='__main__':
    fname = "data/spring_model.nusa"
    m1 = read_spring_model(fname)
    m1.solve()
    m1.simple_report()