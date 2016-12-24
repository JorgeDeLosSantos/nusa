# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************

TRUSS_SIMPLE_REPORT = """
==========================
    NuSA Simple Report
==========================

Model: {model_name}
Number of nodes: {nodes}
Number of elements: {elements}

NODAL DISPLACEMENTS
{nodal_displacements}

NODAL FORCES
{nodal_forces}

ELEMENT FORCES
{element_forces}

ELEMENT STRESSES
{element_stresses}


FINITE-ELEMENT MODEL INFO

NODES:
{nodes_info}

ELEMENTS:
{elements_info}

"""
