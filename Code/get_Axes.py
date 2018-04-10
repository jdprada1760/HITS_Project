#########################################################################
# Calculates the radial profile for the axes at each snapshot
#########################################################################

import matplotlib
matplotlib.use('Agg')
import sys
import gc

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
#import yt
#from yt.units import parsec, Msun
from Py_Libs.shape import *
import os as os

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
#lvl = 'level3_MHD'
#lvl = 'level3_DM'
#lvl = 'level5_Durham' 
#lvl = 'level2'

lvl = sys.argv[1]
halonums = [int(val) for val in sys.argv[2:]]

#halonums = range(1,31)
#halonums = [24]
for j in halonums:
    print("--------------------------------------------------------------------------")
    print(j)
    #halo = "halo" + str(j) +"_MHD"        
    halo = 'halo_'+str(j)
    
    # Considerations of level
    snapnum = 0
    step = 0
    if lvl == 'level3_MHD':
        snapnum = 63
        step = 1
    elif lvl == 'level3_DM':
        snapnum = 127
        step = 2
    
    # This boolean determines the stop of the snapshot calculations
    boole = True
    
    while(boole):
        
        # Loads particles positions, virial radius and other important info
        pos,rvir,info = loadSim(lvl,halo,snapnum)
        pos = np.array(pos,dtype = np.float)
        
        # Stops if 
        if (info['Redshift'] >= 2) : 
            boole = False
        
        print("---------------------------------------------------------------------------------")
        print("----------------------------------------------------------Halo:  "+halo)
        print("----------------------------------------------------------Snap:  "+str(snapnum))
        print("----------------------------------------------------------Rvir:  "+str(rvir))
        print("-------------------------------------------------------NUMPART:  " +str(len(pos)))
        print("---------------------------------------------------------------------------------")
        # Logspace for radius
        logdelta = 2.0/50.0
        #logdelta = 1.0/50.0
        logr_ini = 0.5
        
        # Keeps semiaxes
        semmiaxes = []
        axes_vecs = []
        # Initial radius, axes and inertia
        x_rad = 10.0**(logr_ini)
        axes = np.ones(3,dtype = np.float)*x_rad
        vecs = np.identity(3,dtype = np.float)
        
        for i in range(1000):
        
            print("_______________RAD    : " +str(x_rad))
         
            # Ensures that we obtain an output radius near the corresponding radius
            rad = 0
            l = 0
            while True:
                
                inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))
                rad = (abs(axes[0]*axes[1]*axes[2]))**(1./3.)
                # Advance
                if axes[0] <= 0:
                    break
                elif abs(rad-x_rad)/x_rad > 0.05 :
                    axes = (axes[0]+.5*(x_rad-rad))*axes/axes[0]
                    #axes = np.ones(3,dtype = np.float)*(axes[0]+1.*(x_rad-rad))
                    #vecs = np.identity(3,dtype = np.float)  
                else:
                    break
                    
                l += 1   
                print("+++"+str(l)+"+++",rad,x_rad,axes[0])  
                if l > 30 :
                    break

            x_rad = 10.0**(logr_ini+(i+1)*logdelta)
            if axes[0] > 0:
                semmiaxes.append(np.copy(axes))
                axes_vecs.append(np.copy(vecs).flatten())
                # Uses last result to improve convergence
                axes[1:] = (x_rad/axes[0])*axes[1:]
                axes[0] = x_rad
                #m = 20
                #if i%m == 0:
                    #title = lvl+"_"+halo+"_"+str(i//m)+".png"
                #a,b,c = axes
                    #plotHalo(title,a,b,c,pos,lvl,halo)
            else:
                axes = np.ones(3,dtype = np.float)*x_rad
                vecs = np.identity(3,dtype = np.float)

            #rad = (abs(axes[0]*axes[1]*axes[2]))**(1./3.)
            #print("________________________________________________________________")
            #print(x_rad,rad)
            if( x_rad > 1.0*rvir ):
                break
                
        semmiaxes.append(np.array([rvir,rvir,rvir]))
        semmiaxes = np.array(semmiaxes)
        path = "../Plots/"+lvl+"/"+halo+"/"+"snap_"+str(snapnum)
        # Tries to create directory
        try:
            os.makedirs(path)
        except OSError:
            print("Directory already exists")
        np.savetxt(path+"/s"+str(snapnum)+".txt",semmiaxes, delimiter = ',')
        np.savetxt(path+"/s"+str(snapnum)+"_vecs.txt",axes_vecs, delimiter = ',')
        np.save(path+"/info.npy",info)
        snapnum = snapnum-step
        # Forces garbage collector
        gc.collect()
        

    #np.savetxt("../Plots/"+lvl+"/"+halo+"/rand_sample/"+str(k)+".txt",semmiaxes, delimiter = ',')
    #np.savetxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt",semmiaxes, delimiter = ',')
'''
    if(len(semmiaxes) != 0):
        a,b,c = semmiaxes.T
        yvals = np.array([b/a,c/a,c/b])
        ylabel = ['b/a','c/a','c/b']
        xvals = (a*b*c)**(1./3.)
        fig,axs = ratiosPlots(xvals, yvals, ylabel, rvir)  
        plt.savefig("../Plots/"+lvl+"/"+halo+"/"+"Axial_ratios_vs_R__"+lvl+"_"+halo+".png")
        plt.close()
'''



'''
# Random sampling 10 times
    for k in range(10):
        print("Randomsample___halo:    "+str(k)+"    "+halo)
        pos = np.random.permutation(pos1)
        np.random.shuffle(pos)
        np.random.shuffle(pos)
        pos = pos[:int(len(pos)*1.0/8)]
        print("Radius of the simulation:  "+str(rvir))
        # Logspace for radius
        logdelta = 1.0/50.0
        logr_ini = -1
        # Keeps semiaxes
        semmiaxes = []
        x_rad = 10.0**(logr_ini)
        axes = np.ones(3,dtype = np.float)*x_rad
        vecs = np.identity(3,dtype = np.float)
        for i in range(1000):
            inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))
            x_rad = 10.0**(logr_ini+(i+1)*logdelta)
            if axes[0] > 0:
                semmiaxes.append(np.copy(axes))
                # Uses last result to improve convergence
                axes[1:] = (x_rad/axes[0])*axes[1:]
                axes[0] = x_rad
                #m = 20
                #if i%m == 0:
                    #title = lvl+"_"+halo+"_"+str(i//m)+".png"
                #a,b,c = axes
                    #plotHalo(title,a,b,c,pos,lvl,halo)
            else:
                axes = np.ones(3,dtype = np.float)*x_rad
                vecs = np.identity(3,dtype = np.float)

            rad = (abs(axes[0]*axes[1]*axes[2]))**(1./3.)
            #print("________________________________________________________________")
            #print(x_rad,rad)
            if( x_rad > 2.0*rvir ):
                break
                
'''




