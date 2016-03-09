# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
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
        if "E" in kwargs.keys():
            self.E = kwargs["E"]
    
    def __str__(self):
        return self.name


if __name__=='__main__':
    pass
    
