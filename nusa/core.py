# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************

import numpy as np
import numpy.linalg as la
import _templates as tmp

class Model(object):
	"""
	Model for finite element analysis
	"""
	def __init__(self,name="Model 01"):
		self.name = name
		self.nodes = []
		self.elements = []
		self.dof = 1
		
	def addNode(self,node):
		self.nodes.append(node)
		
	def getNumberOfNodes(self):
		return len(self.nodes)
		
	def getNumberOfElements(self):
		return len(self.elements)
	
	def addElement(self,element):
		"""
		Add element to current model
		
		*element* :  :class:`~core.base.Element`
			Element instance 
		
		::
		
			m1 = Model()
			e1 = Bar(Node((0,0),0),Node((1,0),0))
			m1.addElement(e1)
		
		"""
		self.elements.append(element)
		
	def buildForcesVector(self):
		self.F = np.zeros((self.dof*self.getNumberOfNodes(),1))
		
	def buildGlobalMatrix(self):
		msz = (self.dof)*self.getNumberOfNodes()
		self.KG = np.zeros((msz,msz))
		for element in self.elements:
			ku = element.getElementStiffness()
			n1,n2 = element.getNodes()
			self.KG[n1.label, n1.label] += ku[0,0]
			self.KG[n1.label, n2.label] += ku[0,1]
			self.KG[n2.label, n1.label] += ku[1,0]
			self.KG[n2.label, n2.label] += ku[1,1]
		self.buildForcesVector()
		self.buildDisplacementsVector()
		
	def buildDisplacementsVector(self):
		msz = (self.dof)*self.getNumberOfNodes()
		self.mdl_displacements = dict()
		for label in range(self.getNumberOfNodes()):
			self.mdl_displacements[label] = np.nan
		
	def addForce(self,node,force):
		self.F[node.label] = force[0]
		
	def addConstraint(self,node,**constraint):
		label = node.label
		if constraint.has_key('ux'):
			val = constraint.get('ux')
			node.setDisplacements(ux=val)
			self.mdl_displacements[label] = val
		todel = [lbl for lbl in self.mdl_displacements.keys() \
				if not self.mdl_displacements[lbl] is np.nan]
		self.K2S = np.delete(np.delete(self.KG,todel,0),todel,1)
		self.F2S = np.delete(self.F,todel,0)
		
	def solve(self):
		self.U = la.solve(self.K2S,self.F2S)
		n = 0
		for dsp in self.mdl_displacements.items():
			if self.mdl_displacements[dsp[0]] is np.nan:
				self.mdl_displacements[dsp[0]] = float(self.U[n])
				n += 1
	
	def getElementForces(self):
		for elm in self.elements:
			CK = elm.getElementStiffness()
			n1,n2 = elm.getNodes()
			un1 = self.mdl_displacements[n1.label]
			un2 = self.mdl_displacements[n2.label]
			FV = np.array([un1,un2])
			print np.dot(CK,FV)
		
	def getNodeForces(self):
		return np.dot(self.KG,np.array(self.mdl_displacements.values()))

	def report(self,out="cli"):
		"""
		Write a mini-report for results.
		
		If out is 'cli' then show information in python console.
		"""
		
		if out=="cli":
			print(self)
		elif out=="gui":
			pass
		elif out=="file":
			pass
		
	def __str__(self):
		"""
		String for Model information
		"""
		_str = tmp.MINI_REPORT_TEMPLATE.format(
		model=self.name,
		etype="Truss",
		nelements=self.getNumberOfElements(),
		nnodes=self.getNumberOfNodes(),
		displacements=self.mdl_displacements.values())
		return _str
		
		
class Element(object):
	"""
	Superclass for Element objects 
	"""
	def __init__(self,etype):
		self.etype = etype



class Node(object):
	"""
	Class for node object.
	
	*coordinates* : `tuple`, `list`
		Coordinates of node
	
	*label* : int
		Label of node
		
	::
	
		n1 = Node((0,0),0)
		n2 = Node((0,0),1)
	
	"""
	def __init__(self,coordinates,label):
		self.coordinates = coordinates
		self.__label = label
		self.__ux = np.nan
		self.__uy = np.nan
		
	@property
	def label(self):
		return self.__label
		
	@label.setter
	def label(self,val):
		"""
		Experimental setter for adjust range of labels
		"""
		self.__label = val
		
	@property
	def ux(self):
		return self.__ux
	
	@ux.setter
	def ux(self,val):
		if val in [int,float]:
			self.__ux = val
		else:
			raise ValueError("Value must be float or int")
	
	@property
	def uy(self):
		return self.__uy
	
	@uy.setter
	def uy(self,val):
		if val in [int,float]:
			self.__uy = val
		else:
			raise ValueError("Value must be float or int")
		
	def getLabel(self):
		return self.__label
	
	def getDisplacements(self):
		return (self.ux,self.__uy)
		
	def setDisplacements(self,ux=np.nan,uy=np.nan):
		self.__ux = ux
		self.__uy = uy



class Spring(Element):
	"""
	Bar element for finite element analysis
	
	*nodes* : :class:`~core.base.Node`
		Conectivity for element
	
	*E* : float
		Young modulus
		
	*A* : float
		Area of element
		
	*L* : float
		Length of element
	
	
	Example::
	
		n1 = Node((0,0),1)
		n2 = Node((1,0),1)
		e1 = Spring((n1,n2),200e9,0.02,1)
	
	"""
	def __init__(self,nodes,ke):
		Element.__init__(self,etype="spring")
		self.nodes = nodes
		self.ke = ke
		
	def getElementStiffness(self):
		r"""
		Get stiffness matrix for this element
		
		The stiffness matrix for a spring element is defined by:
		
		.. math::
		
			[K]_e = \left[\begin{matrix}
			k  & -k \\
			-k &  k \\
			\end{matrix} \right]
		
		Return a numpy array.
		"""
		self._KE = np.array([[self.ke,-self.ke],[-self.ke,self.ke]])
		return self._KE
		
	def getNodes(self):
		return self.nodes


class Bar(Element):
	"""
	Bar element for finite element analysis
	
	*nodes* : :class:`~core.base.Node`
		Conectivity for element
	
	*E* : float
		Young modulus
		
	*A* : float
		Area of element
		
	*L* : float
		Length of element
	
	
	Example::
	
		n1 = Node((0,0),1)
		n2 = Node((1,0),1)
		e1 = Bar((n1,n2),200e9,0.02,1)
	
	"""
	def __init__(self,nodes,E,A,L):
		Element.__init__(self,etype="bar")
		self.nodes = nodes
		self.E = E
		self.A = A
		self.L = L
		
	def getElementStiffness(self):
		"""
		Get stiffness matrix for this element
		
		"""
		self._K = (self.A*self.E/self.L)*np.array([[1,-1],[-1,1]])
		return self._K
		
	def getNodes(self):
		return self.nodes
		

if __name__=='__main__':
	pass
