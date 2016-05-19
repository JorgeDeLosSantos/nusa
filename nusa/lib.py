# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
"""
Defining some common materials.
"""

import material as mat

## ================== Materials ======================
STEEL_1018 = mat.Material("1018 Steel",E=205e9,nu=0.3)
ALUMINIUM_6061 = mat.Material("6061 Aluminium Alloy",E=69e9,nu=0.33)
