# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.blogspot.mx
#  License: MIT License
# ***********************************

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


if __name__=='__main__':
    a = Material("Aluminio",E=200,nu=0.3)
    print dir(a)
    
