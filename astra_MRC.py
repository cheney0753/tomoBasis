#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 16:45:01 2016

@author: zhong
"""
import numpy as np
import astra
import struct

def read(fname, returnId= False):

#Read the file

    fn= open(fname)
    
    dims= struct.unpack("3I", fn.read(3*4))

    dType= struct.unpack("I", fn.read(4))[0]
    
    fn.seek( 23*4, 0)
    
    nextId=  struct.unpack("I", fn.read(4))[0]
    
    angles= np.zeros( dims[2])
    
    
    for i in range(dims[2]):
        fn.seek(1024+i*128, 0)
        angles[i] = struct.unpack("f", fn.read(4))[0]
        
    angles= angles/180*np.pi
    

    proj_geom = astra.create_proj_geom('parallel3d', 1,1, dims[1], dims[0], angles)
    
    if dType == 1:
        dataType= 'h' #the datatype is short (int32 in MATLAB)
    elif dType == 2:
        dataType= 'f' #the datatype is float (single in MATLAB)
    else:
        print("Data type error: only support short integer or float")
        raise Exception
        
    V= np.zeros((dims[1], dims[2], dims[0]) ,dtype=np.float32)
    
    fn.seek(1024+nextId, 0)
    
    for i in range(dims[2]):
        
        if dType == 1:
            V[:,i,:]= np.asarray(struct.unpack( str(dims[0]*dims[1])+dataType, \
                        fn.read(dims[0]*dims[1]*2))).reshape((dims[0], dims[1]))
        elif dType == 2:
            V[:,i,:]= np.asarray(struct.unpack( str(dims[0]*dims[1])+dataType, \
                        fn.read(dims[0]*dims[1]*4))).reshape((dims[0], dims[1]))
        else:
            print("Data type error: only support short integer or float")
            raise Exception
            
    V -= np.min(V)
            
    proj_id= astra.data3d.create('-sino', proj_geom, V)

    
    # show the images
    import pylab
    
    pylab.figure(1)
    pylab.imshow(V[:,10,:])
    pylab.show()
    
    if returnId:
        return V, proj_geom, proj_id
    else:
        astra.data3d.delete(proj_id)
        return V, proj_geom

def write(fname, data, geometry, precision='float'):
    
#Define the precision: only support float and short
   
    if precision == 'short':
        dType= 1
    elif precision == 'float':
        dType= 2
    else:
        print("Data type error: only support short or float")
        raise Exception
        
#Open the file
    try:
        fn= open(fname,  'wb')
    except:
        raise Exception("Could not open "+ fname +" for writing.")

#Define the dimension
    if geometry.has_key('GridSliceCount'):
        geomType= 'Vol'
    elif geometry.has_key('DetectorRowCount'):
        geomType= 'Proj3D'
    elif geometry.has_key('DetectorCount'):
        geomType= 'Proj2D'
    
        
#Define the geometry
    if geomType == 'Vol':
        dims=geometry['GridColCount']*geometry['GridRowCount']* geometry['GridSliceCount']
    elif geomType == 'Proj3D':
        dims=[geometry['DetectorColCount'], geometry['DetectorRowCount'],  len(geometry['ProjectionAngles'])]
    elif geomType == 'Proj2D':
        dims=[1,  geometry['DetectorCount'], len(geometry['ProjectionAngles'])]
    
#Define the file size
    if dType == 1:  
        sizePack= 2
    else: 
        sizePack= 4
        
    byteN= np.prod(dims)*sizePack
    byteN+=1024 #Add size of the header
    
    #Add extended header
    if geomType== 'Proj2D' or geomType == 'Proj3D':
        byteN += 128 * len(geometry['ProjectionAngles'])
          
    #Pre-allocate (??)
    fn.write(struct.pack( str(byteN) +'b', *(np.zeros(byteN, dtype=np.int8).tolist() ) ))
    
    #Go back to the beginning
    fn.seek(0, 0)
    
    #Write dimensions
    fn.write(struct.pack( str(len(dims)) +'I', *dims))
    
    #Write mode
    if dType == 1:
        fn.write(struct.pack('I', 1))
    elif dType == 2:
        fn.write(struct.pack('I', 2))
    
    #store grid size (??)
    fn.seek(28, 0)
    fn.write(struct.pack(str(len(dims)) +'I', *dims))
    
    #Get cell size
    if geomType== 'Proj2D':
        nx= geometry['DetectorWidth']     
        ny= 1
    elif geomType == 'Proj3D':
        nx= geometry['DetectorSpacingX'] 
        ny= geometry['DetectorSpacingY']
    elif geomType == 'Vol':
        nx= 1
        ny= 1
    
    nz= 1

    #Write cell size
    fn.write(struct.pack( '3I', nx, ny, nz))
    fn.seek(64, 0)
    
    #Need for IMOD
    fn.write(struct.pack('3I', 1, 2, 3))
    
    #For scaling
    amin= np.min(data)    
    amax= np.max(data)
    amean= np.mean(data)
    
    fn.write(struct.pack('3f', amin, amax, amean))
    
    #Write the angles
    if geomType== 'Proj2D' or geomType == 'Proj3D':
        
        nAngles= len(geometry['ProjectionAngles'])
        
        #Skip the bytes in the header
        fn.seek(92,0)
        nextId= nAngles * 128
        fn.write(struct.pack('I', nextId) )
        
        angles= geometry['ProjectionAngles']*180/np.pi
        
        #Write the angles
        for i in range(len(angles)):
            fn.seek(1024+i*128, 0)
            fn.write(struct.pack('f', angles[i]))
            
    else: #No angles 
        fn.seek(92, 0) #Extended header is empty
        fn.write(struct.pack('I', 0))
        nextId= 0;

#Now we write the actual data
    fn.seek(1024+nextId, 0)
    
    if dType == 1:
        data-= 32768
        
    #Permute the data if it is projection data
    if geomType == 'Proj2D' or geomType == 'Proj3D':
        data=np.transpose(data, (0, 2, 1))
        
    if dType == 1:
        dataType= 'I'
    else:
        dataType= 'f'
    
    for i in range(data.shape[2]):
        fn.write(struct.pack(str( data.shape[0]*data.shape[1]) + dataType,*( data[:,:,i].reshape(data.shape[0]*data.shape[1],1) ) ))
        
    return