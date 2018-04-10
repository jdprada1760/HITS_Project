#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import gc

#########################################################################
# Calculates the axes in terms of the virial radius at each redshift
#########################################################################

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import yt
from yt.units import parsec, Msun
from Py_Libs.shape import *
import os as os
import sys

# C imports
import ctypes as ct
inertia = ct.cdll.LoadLibrary('./C_Libs/inertia.so')

# Simiulation specs
#lvl = 'level5'
#lvl = 'level4_MHD'
#lvl = 'level5_Durham' 
#halonums = [27]
#halonums = range(1,11)
#halonums = ['halo16','halo16_MHD','halo24','halo24_MHD','halo28', 'halo6_MHD','halo9','halo9_MHD']

'''
# Loads simulation
snapnum = 127
pos,rvir = loadSim(lvl,'halo_1',snapnum)
pos = np.array(pos,dtype = np.float)
eigenvals = rvir*np.ones(3,dtype = np.float)
eigenvecs = np.identity(3,dtype = np.float)
inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(eigenvecs.ctypes.data), ct.c_void_p(eigenvals.ctypes.data))
print(eigenvals)
print(eigenvecs)
'''

#lvl = 'level5'
#lvl = 'level3_MHD' #
#lvl = 'level3_DM' #halo21
#lvl = 'level5_Durham' 
#lvl = 'level2'
lvl = sys.argv[1]
halonums = [6,16,21,23,24,27]
#halonums = range(1,31)
#halonums = [24]
for j in halonums:
    print("--------------------------------------------------------------------------")
    print(j)
    #halo = "halo" + str(j) +"_MHD"        
    halo = 'halo_'+str(j)
    
    #snapnum = 17
    #snapnum = 127
    #snapnum = 63
    snapnum = int(sys.argv[2])    
    # This boolean determines the stop of the snapshot calculations
    boole = True
    
    # Saves all the information needed
    semmiaxes = []
    axes_vecs = []
    
    pos,rvirrr,info = loadSim(lvl,halo,63)
    while(boole):
        
        # Loads particles positions, virial radius and other important info
        pos,rvir,info = loadSim(lvl,halo,snapnum)
        # We need the virial radius in physical coordinates
        rvir = rvirrr*info["Scale factor"]
        pos = np.array(pos,dtype = np.float)
        
        # Stops if 
        if (rvir < 10) : 
            boole = False
        
        #print("---------------------------------------------------------------------------------")
        #print("----------------------------------------------------------Halo:  "+halo)
        #print("----------------------------------------------------------Snap:  "+str(snapnum))
        #print("----------------------------------------------------------Rvir:  "+str(rvir))
        #print("-------------------------------------------------------NUMPART:  " +str(len(pos)))
        #print("---------------------------------------------------------------------------------")
       
        #print("_______________RAD    : " +str(rvir))
     
        axes = np.ones(3,dtype = np.float)*rvir
        vecs = np.identity(3,dtype = np.float)

        # Ensures that we obtain an output radius near the corresponding radius
        rad = 0
        l = 0
        while True:
            
            inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))
            rad = (abs(axes[0]*axes[1]*axes[2]))**(1./3.)
            # Advance
            if axes[0] < 0:
                break
            elif abs(rad-rvir)/rvir > 0.01 :
                axes = (axes[0]+0.5*(rvir-rad))*axes/axes[0]
		axes_vecs.append(np.copy(vecs).flatten())
            else:
                break
                
            l += 1   
            #print("+++it+++,rad,rvir,axes[0]")
            #print("+++"+str(l)+"+++",rad,rvir,axes[0])  
            if l > 200 :
                break

        if axes[0] > 0:
            axx = np.copy(axes)
            semmiaxes.append([info['Redshift'],axx[0],axx[1],axx[2]])
            
        snapnum = snapnum-1 
        gc.collect()
    semmiaxes = np.array(semmiaxes)
    path = "../Plots/"+lvl+"/"+halo
    np.savetxt(path+"/Z_axes"+".txt",semmiaxes, delimiter = ',')
    np.savetxt(path+"/Z_vecs"+".txt",axes_vecs, delimiter = ',')       

    #np.savetxt("../Plots/"+lvl+"/"+halo+"/rand_sample/"+str(k)+".txt",semmiaxes, delimiter = ',')
    #np.savetxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt",semmiaxes, delimiter = ',')

