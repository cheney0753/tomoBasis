#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 16:45:01 2016

@author: zhong
"""

def astra_read_FEI_MRC(fname):
    import numpy as np
    import astra
    import struct
#Read the file

    fn = open(fname)
    
    dims= struct.unpack("3I", fn.read(3*4))

    dType = struct.unpack("I", fn.read(4))[0]
    
    fn.seek( 23*4, 0)
    
    nextId =  struct.unpack("I", fn.read(4))[0]
    
    angles = np.zeros( dims[2])
    
    
    for i in range(dims[2]):
        fn.seek(1024+i*128, 0)
        angles[i] = struct.unpack("f", fn.read(4))[0]
        
    angles= angles/180*np.pi
    

    proj_geom = astra.create_proj_geom('parallel3d', 1,1, dims[1], dims[0], angles)
    
    if dType == 1:
        dataType = 'h' #the datatype is short (int32 in MATLAB)
    elif dType == 2:
        dataType = 'f' #the datatype is float (sing;e in MATLAB)
    else:
        print "Data type error: only support short integer or float"
        raise Exception
        
    V= np.zeros((dims[1], dims[2], dims[0]) ,dtype=np.float32)
    
    fn.seek(1024+nextId, 0)
    
    for i in range(dims[2]):
        
        if dType == 1:
            V[:,i,:] = np.asarray(struct.unpack( str(dims[0]*dims[1])+dataType, \
                        fn.read(dims[0]*dims[1]*2))).reshape((dims[0], dims[1]))
        elif dType == 2:
            V[:,i,:] = np.asarray(struct.unpack( str(dims[0]*dims[1])+dataType, \
                        fn.read(dims[0]*dims[1]*4))).reshape((dims[0], dims[1]))
        else:
            print "Data type error: only support short integer or float"
            raise Exception
            
    V -= np.min(V)
            
    proj_id= astra.data3d.create('-sino', proj_geom, V)

    
    # show the images
    #import pylab
    
    #pylab.figure(1)
    #pylab.imshow(V[:,10,:])
    #pylab.show()
    
    return [proj_id, V]
