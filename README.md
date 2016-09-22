# NuSA

Numerical Structural Analysis with Python

## Version / Status

Current: **0.1.0-dev**

## Requirements

* NumPy
* Matplotlib (Optional for post)


## Installation

Using pip:

```
$ pip install git+https://github.com/JorgeDeLosSantos/nusa.git
```

## Capabilities

* Bar elements analysis
* 2D frames analysis
* Beams analysis

## Mini-Demo

**Example 01**. For the spring assemblage with arbitrarily numbered nodes shown in the figure 
obtain (a) the global stiffness matrix, (b) the displacements of nodes 3 and 4, (c) the 
reaction forces at nodes 1 and 2, and (d) the forces in each spring. A force of 5000 lb
is applied at node 4 in the $x$ direction. The spring constants are given in the figure.
Nodes 1 and 2 are fixed.

![](docs\nusa-theory\src\spring-element\example_01.PNG)

```python
# -*- coding: utf-8 -*-
# NuSA Demo
from nusa.core import *
from nusa.model import *
from nusa.element import *
    
def test1():
    """
    Logan, D. (2007). A first course in the finite element analysis.
    Example 2.1, pp. 42.
    """
    P = 5000.0

    # Model
    m1 = SpringModel("2D Model")

    # Nodes
    n1 = Node((0,0))
    n2 = Node((0,0))
    n3 = Node((0,0))
    n4 = Node((0,0))

    # Elements
    e1 = Spring((n1,n3),1000.0)
    e2 = Spring((n3,n4),2000.0)
    e3 = Spring((n4,n2),3000.0)

    # Add elements 
    for nd in (n1,n2,n3,n4):
        m1.addNode(nd)
    for el in (e1,e2,e3):
        m1.addElement(el)

    m1.buildGlobalMatrix()
    m1.addForce(n4, (P,))
    m1.addConstraint(n1, ux=0)
    m1.addConstraint(n2, ux=0)
    m1.solve()

if __name__ == '__main__':
    test1()
```



## Documentation



## About...

```
Developer: Pedro Jorge De Los Santos
E-mail: delossantosmfq@gmail.com
Blog: labdls.blogspot.mx
```