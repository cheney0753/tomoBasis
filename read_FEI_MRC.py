#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 16:45:01 2016

@author: zhong
"""

def astra_read_FEI_MRC(fname):
    
    import numpy as np
#Read the file
    fn = read(fname)
    
    dims = np.zeros(3,1, int)
    for i in range(0,2):
        dims(i) = struct.unpack("I", fn.read(4))

    dimType = struct.unpack("I", fn.read(4))
    
    seek(fn, [23*4, 0])
    
    nextId =  struct.unpack("I", fn.read(4))
    
    