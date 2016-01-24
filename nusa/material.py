# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************

class Material(object):
	def __init__(self,name,**kwargs):
		self.name = name
		for arg in kwargs.keys():
			if arg is "E":
				self.E = kwargs[arg]
			else:
				pass
	
	def getE(self):
		pass
		


if __name__=='__main__':
	m1 = Material("Acero",E=200e9)
	print m1.E
	
