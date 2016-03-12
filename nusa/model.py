# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import numpy as np
import numpy.linalg as la
import templates as tmp
from core import Model
from _info import __version__,__organization__

class SpringModel(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="Spring Model 01"):
        Model.__init__(self,name=name,mtype="spring")
        self.F = {}
        self.U = {}
        self.dof = 1 # 1 DOF for spring element
        
    def buildForcesVector(self):
        for node in self.nodes.values():
            self.F[node.label] = 0
        
    def buildGlobalMatrix(self):
        msz = (self.dof)*self.getNumberOfNodes()
        self.KG = np.zeros((msz,msz))
        for element in self.elements.values():
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
        for node in self.nodes.values():
            self.U[node.label] = np.nan
        
    def addForce(self,node,force):
        self.F[node.label] = force[0]
        
    def addConstraint(self,node,**constraint):
        label = node.label
        if constraint.has_key('ux'):
            val = constraint.get('ux')
            node.setDisplacements(ux=val)
            self.U[label] = val
        
    def solve(self):
        knw = [lbl for lbl in self.U.keys() if not self.U[lbl] is np.nan]
        unknw = [lbl for lbl in self.U.keys() if self.U[lbl] is np.nan]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.F.values(),knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        self.U.update(zip(unknw,self.solved_u))
        # For nodal forces/reactions
        self.NF = self.F.copy()
        nf_calc = np.dot(self.KG,self.U.values())
        self.NF.update(zip(self.NF.keys(),nf_calc))
        # For element forces
        self.EF = self.elements.copy()
        for elm in self.elements.values():
            _ke = elm.getElementStiffness()
            na,nb = elm.getNodes()
            una = self.U[na.label]
            unb = self.U[nb.label]
            _uv = np.array([una,unb])
            self.EF.update([(elm.label,np.dot(_ke,_uv)),])
            
    def getDisplacements(self):
        _str = ""
        for lbl,u in self.U.items():
            _str += (str(lbl+1) + "\t|\t" + str(u) + "\n")
        return _str
                
    def getNodalForces(self):
        _str = ""
        for lbl,nf in self.NF.items():
            _str += (str(lbl+1) + "\t|\t" + str(nf) + "\n")
        return _str
    
    def getElementForces(self):
        _str = ""
        for lbl,ef in self.EF.items():
            _str += (str(lbl+1) + "\t|\t" + str(ef) + "\n")
        return _str
        
    def report(self,out="cli"):
        """
        Write a mini-report for results.
        
        If out is 'cli' then show information in python console.
        """
        if out=="cli":
            print(self.CLIReport())
        elif out=="gui":
            pass
        elif out=="file":
            pass
        
    def CLIReport(self):
        """
        String for Spring Model information
        """
        _str = tmp.MINI_REPORT_TEMPLATE.format(
        model=self.name,
        mtype=self.mtype,
        nelements=self.getNumberOfElements(),
        nnodes=self.getNumberOfNodes(),
        displacements=self.getDisplacements(),
        nodalforces=self.getNodalForces(),
        elementforces=self.getElementForces(),
        version=__version__)
        return _str
        

class BarModel(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="Bar Model 01"):
        Model.__init__(self,name=name,mtype="bar")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 1 # 1 DOF for bar element
        
    def buildForcesVector(self):
        for node in self.nodes.values():
            self.F[node.label] = 0
        
    def buildGlobalMatrix(self):
        msz = (self.dof)*self.getNumberOfNodes()
        self.KG = np.zeros((msz,msz))
        for element in self.elements.values():
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
        for node in self.nodes.values():
            self.U[node.label] = np.nan
        
    def addForce(self,node,force):
        self.F[node.label] = force[0]
        
    def addConstraint(self,node,**constraint):
        label = node.label
        if constraint.has_key('ux'):
            val = constraint.get('ux')
            node.setDisplacements(ux=val)
            self.U[label] = val
        
    def solve(self):
        knw = [lbl for lbl in self.U.keys() if not self.U[lbl] is np.nan]
        unknw = [lbl for lbl in self.U.keys() if self.U[lbl] is np.nan]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.F.values(),knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        self.U.update(zip(unknw,self.solved_u))
        # For nodal forces/reactions
        self.NF = self.F.copy()
        nf_calc = np.dot(self.KG,self.U.values())
        self.NF.update(zip(self.NF.keys(),nf_calc))
        # For element forces
        self.EF = self.elements.copy()
        for elm in self.elements.values():
            _ke = elm.getElementStiffness()
            na,nb = elm.getNodes()
            una = self.U[na.label]
            unb = self.U[nb.label]
            _uv = np.array([una,unb])
            self.EF.update([(elm.label,np.dot(_ke,_uv)),])
            
    def getDisplacements(self):
        _str = ""
        for lbl,u in self.U.items():
            _str += (str(lbl+1) + "\t|\t" + str(u) + "\n")
        return _str
                
    def getNodalForces(self):
        _str = ""
        for lbl,nf in self.NF.items():
            _str += (str(lbl+1) + "\t|\t" + str(nf) + "\n")
        return _str
    
    def getElementForces(self):
        _str = ""
        for lbl,ef in self.EF.items():
            _str += (str(lbl+1) + "\t|\t" + str(ef) + "\n")
        return _str
        
    def report(self,out="cli",filename="report.txt"):
        """
        Write a mini-report for results.
        
        If out is 'cli' then show information in python console.
        """
        if out=="cli":
            print(self.CLIReport())
        elif out=="gui":
            pass
        elif out=="file":
            self.fileReport(filename)
        
    def CLIReport(self):
        """
        String for Spring Model information
        """
        _str = tmp.MINI_REPORT_TEMPLATE.format(
        model=self.name,
        mtype=self.mtype,
        nelements=self.getNumberOfElements(),
        nnodes=self.getNumberOfNodes(),
        displacements=self.getDisplacements(),
        nodalforces=self.getNodalForces(),
        elementforces=self.getElementForces(),
        version=__version__)
        return _str
        
    def fileReport(self,filename):
        """
        Write report to file
        """
        out_str = self.CLIReport()
        f = open(filename,"w")
        f.write(out_str)
        f.close()


class BeamModel(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="Beam Model 01"):
        Model.__init__(self,name=name,mtype="beam")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 2 # 2 DOF for beam element
        
    def buildForcesVector(self):
        for node in self.nodes.values():
            self.F[node.label] = 0
        
    def buildGlobalMatrix(self):
        msz = (self.dof)*self.getNumberOfNodes()
        self.KG = np.zeros((msz,msz))
        for element in self.elements.values():
            ku = element.getElementStiffness()
            n1,n2 = element.getNodes()
            self.KG[2*n1.label, 2*n1.label] += ku[0,0]
            self.KG[2*n1.label, 2*n1.label+1] += ku[0,1]
            self.KG[2*n1.label, 2*n2.label] += ku[0,2]
            self.KG[2*n1.label, 2*n2.label+1] += ku[0,3]
            
            self.KG[2*n1.label+1, 2*n1.label] += ku[1,0]
            self.KG[2*n1.label+1, 2*n1.label+1] += ku[1,1]
            self.KG[2*n1.label+1, 2*n2.label] += ku[1,2]
            self.KG[2*n1.label+1, 2*n2.label+1] += ku[1,3]
            
            self.KG[2*n2.label, 2*n1.label] += ku[2,0]
            self.KG[2*n2.label, 2*n1.label+1] += ku[2,1]
            self.KG[2*n2.label, 2*n2.label] += ku[2,2]
            self.KG[2*n2.label, 2*n2.label+1] += ku[2,3]
            
            self.KG[2*n2.label+1, 2*n1.label] += ku[3,0]
            self.KG[2*n2.label+1, 2*n1.label+1] += ku[3,1]
            self.KG[2*n2.label+1, 2*n2.label] += ku[3,2]
            self.KG[2*n2.label+1, 2*n2.label+1] += ku[3,3]
            #~ print element.label, 2*n1.label, 2*n1.label+1, 2*n2.label, 2*n2.label+1
            
        self.buildForcesVector()
        self.buildDisplacementsVector()
        
        
    def buildDisplacementsVector(self):
        msz = (self.dof)*self.getNumberOfNodes()
        for node in self.nodes.values():
            self.U[node.label] = np.nan
        
    def addForce(self,node,force):
        self.F[node.label] = force[0]
        
    def addConstraint(self,node,**constraint):
        label = node.label
        if constraint.has_key('ux'):
            val = constraint.get('ux')
            node.setDisplacements(ux=val)
            self.U[label] = val
        if constraint.has_key('uy'):
            val = constraint.get('ux')
            node.setDisplacements(uy=val)
            self.U[label] = val
        
    def solve(self):
        knw = [lbl for lbl in self.U.keys() if not self.U[lbl] is np.nan]
        unknw = [lbl for lbl in self.U.keys() if self.U[lbl] is np.nan]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.F.values(),knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        self.U.update(zip(unknw,self.solved_u))
        # For nodal forces/reactions
        self.NF = self.F.copy()
        nf_calc = np.dot(self.KG,self.U.values())
        self.NF.update(zip(self.NF.keys(),nf_calc))
        # For element forces
        self.EF = self.elements.copy()
        for elm in self.elements.values():
            _ke = elm.getElementStiffness()
            na,nb = elm.getNodes()
            una = self.U[na.label]
            unb = self.U[nb.label]
            _uv = np.array([una,unb])
            self.EF.update([(elm.label,np.dot(_ke,_uv)),])
            
    def getDisplacements(self):
        _str = ""
        for lbl,u in self.U.items():
            _str += (str(lbl+1) + "\t|\t" + str(u) + "\n")
        return _str
                
    def getNodalForces(self):
        _str = ""
        for lbl,nf in self.NF.items():
            _str += (str(lbl+1) + "\t|\t" + str(nf) + "\n")
        return _str
    
    def getElementForces(self):
        _str = ""
        for lbl,ef in self.EF.items():
            _str += (str(lbl+1) + "\t|\t" + str(ef) + "\n")
        return _str
        
    def report(self,out="cli",filename="report.txt"):
        """
        Write a mini-report for results.
        
        If out is 'cli' then show information in python console.
        """
        if out=="cli":
            print(self.CLIReport())
        elif out=="gui":
            pass
        elif out=="file":
            self.fileReport(filename)
        
    def CLIReport(self):
        """
        String for Spring Model information
        """
        _str = tmp.MINI_REPORT_TEMPLATE.format(
        model=self.name,
        mtype=self.mtype,
        nelements=self.getNumberOfElements(),
        nnodes=self.getNumberOfNodes(),
        displacements=self.getDisplacements(),
        nodalforces=self.getNodalForces(),
        elementforces=self.getElementForces(),
        version=__version__)
        return _str
        
    def fileReport(self,filename):
        """
        Write report to file
        """
        out_str = self.CLIReport()
        f = open(filename,"w")
        f.write(out_str)
        f.close()



if __name__=='__main__':
    pass
