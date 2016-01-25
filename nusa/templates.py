# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
MINI_REPORT_TEMPLATE = """
==========================
	NuSA {version}    
==========================


Mini-Report for {model}

Model type: {mtype}
Number of nodes: {nnodes}
Number of elements: {nelements}


=========
Solutions
=========

-------------------------------------------------------------
Displacements:

{displacements}


-------------------------------------------------------------
Nodal forces:

{nodalforces}


-------------------------------------------------------------
Element forces:

{elementforces}

"""
