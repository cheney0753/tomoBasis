#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 14:12:40 2016

@author: zhong
"""
# This file is a wrapper for SIRT algorithm in astra-toolbox


import numpy as np
import astra

def SIRT_astra(proj_geom, vol_geom, P, itN, nonnegativity= False ):
    
# The nonnegativity constraint is off in default
# Determine whether 2D or 3D

    if proj_geom['type']=='parallel':
# 2D geometry 
        P= np.reshape(P, len(proj_geom['ProjectionAngles']), proj_geom['DetectorCount']);
    
        proj_id= astra.create_projector('cuda', proj_geom, vol_geom)
        
        sino_id= astra.create_sino (P, proj_id)
        
        rec_id = astra.data2d.create('-vol', vol_geom)
        
        # define the algorithm
        cfg=  astra.astra_dict('SIRT_CUDA');
        cfg['ReconstructionDataId']= rec_id
        cfg['ProjectionDataId']=sino_id

        # set the nonnegativity constraint
        if nonnegativity==True:
            cfg['option']= {}
            cfg['option']['MinConstraint']=0

        alg_id= astra.algorithm.create(cfg)
        astra.algorithm.run(alg_id, itN);       
        
        #retrieve the reconstruction volume
        rec= astra.data2d.get(rec_id)
        
        #delete the data
        astra.data2d.delete(sino_id, rec_id)
        astra.projector.delete(proj_id)
        astra.algorithm.delete(alg_id)
        
        return rec
        
    elif proj_geom['type']=='parallel3d':
# 3D geometry
        P=np.reshape(P, proj_geom['DetectorColCount'], len(proj_geom['ProjectionAngles0']),     proj_geom['DetectorRowCount'])
        
        
        
        
        
        
    else:
        assert('Only parallel and parallel3d are supported!')
        
        
