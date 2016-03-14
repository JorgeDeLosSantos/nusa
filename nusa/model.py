# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import re
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
            self.F[node.label] = {"fy":0, "m":0} # (fy, m)
        
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
        for node in self.nodes.values():
            self.U[node.label] = {"uy":np.nan, "r":np.nan} # (uy, r)
    
    def addForce(self,node,force):
        self.F[node.label]["fy"] = force[0]
        
    def addMoment(self,node,moment):
        self.F[node.label]["m"] = force[0]
        
    def addConstraint(self,node,**constraint):
        cs = constraint
        if cs.has_key('ux') and cs.has_key("uy") and cs.has_key('r'): # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            r = cs.get('r')
            node.setDisplacements(ux=ux, uy=uy)
            print("Empotrado")
            self.U[node.label]["uy"] = uy
            self.U[node.label]["r"] = r
        elif cs.has_key('ux') and cs.has_key("uy"): # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.setDisplacements(ux=ux, uy=uy)
            print("Fijado")
            self.U[node.label]["uy"] = uy
        elif cs.has_key('uy'):
            uy = cs.get('uy')
            node.setDisplacements(uy=uy)
            print("Soporte simple")
            self.U[node.label]["uy"] = uy
        
    def solve(self):
        self.VU = [node[key] for node in self.U.values() for key in ("uy","r")]
        self.VF = [node[key] for node in self.F.values() for key in ("fy","m")]
        knw = [pos for pos,value in enumerate(self.VU) if not value is np.nan]
        unknw = [pos for pos,value in enumerate(self.VU) if value is np.nan]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.VF,knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        for k,ic in enumerate(unknw):
            nd, var = self.index2key(ic)
            self.U[nd][var] = self.solved_u[k]
        # For nodal forces/reactions
        self.NF = self.F.copy()
        self.VU = [node[key] for node in self.U.values() for key in ("uy","r")]
        nf_calc = np.dot(self.KG, self.VU)
        for k,ic in enumerate(range(2*self.getNumberOfNodes())):
            nd, var = self.index2key(ic, ("fy","m"))
            self.NF[nd][var] = nf_calc[k]
        
            
    def index2key(self,idx,opts=("uy","r")):
        node = idx/2
        var = opts[0] if ((-1)**idx)==1 else opts[1]
        return node,var
            
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
        

class TrussModel(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="Truss Model 01"):
        Model.__init__(self,name=name,mtype="truss")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 2 # 2 DOF for truss element
        
    def buildForcesVector(self):
        for node in self.nodes.values():
            self.F[node.label] = {"fx":0, "fy":0} # (fy, m)
        
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
        for node in self.nodes.values():
            self.U[node.label] = {"ux":np.nan, "uy":np.nan} # (uy, r)
    
    def addForce(self,node,force):
        self.F[node.label]["fx"] = force[0]
        self.F[node.label]["fy"] = force[1]
        
    def addConstraint(self,node,**constraint):
        cs = constraint
        if cs.has_key('ux') and cs.has_key("uy"):
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.setDisplacements(ux=ux, uy=uy)
            print("Fijado")
            self.U[node.label]["ux"] = ux
            self.U[node.label]["uy"] = uy
        elif cs.has_key('ux'):
            ux = cs.get('ux')
            node.setDisplacements(ux=ux)
            print("Soporte simple x")
            self.U[node.label]["ux"] = ux
        elif cs.has_key('uy'):
            uy = cs.get('uy')
            node.setDisplacements(uy=uy)
            print("Soporte simple y")
            self.U[node.label]["uy"] = uy
        
    def solve(self):
        self.VU = [node[key] for node in self.U.values() for key in ("ux","uy")]
        self.VF = [node[key] for node in self.F.values() for key in ("fx","fy")]
        knw = [pos for pos,value in enumerate(self.VU) if not value is np.nan]
        unknw = [pos for pos,value in enumerate(self.VU) if value is np.nan]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.VF,knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        for k,ic in enumerate(unknw):
            nd, var = self.index2key(ic)
            self.U[nd][var] = self.solved_u[k]
        # For nodal forces/reactions
        self.NF = self.F.copy()
        self.VU = [node[key] for node in self.U.values() for key in ("ux","uy")]
        nf_calc = np.dot(self.KG, self.VU)
        for k,ic in enumerate(range(2*self.getNumberOfNodes())):
            nd, var = self.index2key(ic, ("fx","fy"))
            self.NF[nd][var] = nf_calc[k]
        
            
    def index2key(self,idx,opts=("ux","uy")):
        node = idx/2
        var = opts[0] if ((-1)**idx)==1 else opts[1]
        return node,var
            
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
        elementforces="",
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


class TrussModel2(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="Truss Model 01"):
        Model.__init__(self,name=name,mtype="truss")
        self.F = [] # Forces
        self.U = [] # Displacements
        self.dof = 2 # 2 DOF for truss element
        
    def buildForcesVector(self):
        for node in self.nodes.values():
            fx, fy = node.getForces()
            self.F.append(fx)
            self.F.append(fy)
        self.F = np.array(self.F)
        
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
            
        #~ self.buildForcesVector()
        #~ self.buildDisplacementsVector()
        
    def buildDisplacementsVector(self):
        for node in self.nodes.values():
            ux, uy = node.getDisplacements()
            self.U.append(ux)
            self.U.append(uy)
        self.U = np.array(self.U)
    
    def addForce(self,node,force):
        node.setForces(fx=force[0], fy=force[1])
        
    def addConstraint(self,node,**constraint):
        cs = constraint
        if cs.has_key('ux') and cs.has_key("uy"):
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.setDisplacements(ux=ux, uy=uy)
        elif cs.has_key('ux'):
            ux = cs.get('ux')
            node.setDisplacements(ux=ux)
        elif cs.has_key('uy'):
            uy = cs.get('uy')
            node.setDisplacements(uy=uy)
        
    def solve(self):
        self.buildDisplacementsVector()
        self.buildForcesVector()
        knw = [pos for pos,value in enumerate(self.U) if not np.isnan(value)]
        unknw = [pos for pos,value in enumerate(self.U) if np.isnan(value)]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.F,knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        for k,idx in enumerate(unknw):
            nd = idx/2
            self.U[idx] = self.solved_u[k]
            if (-1)**idx > 0:
                self.nodes[nd].ux = self.solved_u[k]
            else:
                self.nodes[nd].uy = self.solved_u[k]
        #~ # For nodal forces/reactions
        self.NF = self.F.copy()
        nf_calc = np.dot(self.KG, self.U)
        for nd in self.nodes.keys():
            self.NF[2*nd] = nf_calc[2*nd]
            self.NF[2*nd+1] = nf_calc[2*nd+1]
        print self.NF
            
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
        elementforces="",
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
