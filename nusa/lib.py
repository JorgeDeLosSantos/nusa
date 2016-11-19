# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
"""
Defining some common materials.
"""

class Material(object):
    def __init__(self,name,**kwargs):
        """
        
        Parameters
        ----------
        
        name :    Name of material
        
        
        **kwargs:
        
        E       :    Elastic modulus
        nu      :    Poisson ratio
        density :    Density
        
        """
        self.name = name
        self.__dict__ = kwargs # add all properties
    
    def __str__(self):
        return self.name

## ================== Materials ======================

STEEL_1018 = Material("1018 Steel", E=205e9, nu=0.3)
ALUMINIUM_6061 = Material("6061 Aluminium Alloy", E=69e9, nu=0.33)

