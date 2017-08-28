import matplotlib
matplotlib.use('Agg')

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import yt
from yt.units import parsec, Msun
from Py_Libs.shape import *

# C imports
import ctypes as ct
inertia = ct.cdll.LoadLibrary('./C_Libs/inertia.so')

# Simiulation specs
#lvl = 'level5'
lvl = 'level4_MHD'
#lvl = 'level5_Durham' 
halonums = [20]
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
lvl = 'level4_MHD'
#lvl = 'level5_Durham' 
#halonums = [20]
halonums = range(1,31)
#halonums = ['halo16','halo16_MHD','halo24','halo24_MHD','halo28', 'halo6_MHD','halo9','halo9_MHD']
for j in halonums:
    print("--------------------------------------------------------------------------")
    print(j)
    halo = 'halo_'+str(j)
    #halo = j
    #snap = "/snapdir_127/"
    #subSnap = "/groups_127/"
    snapnum = 127
    #snapnum = 63
    #snapnum = 255
    pos,rvir = loadSim(lvl,halo,snapnum)
    pos = np.array(pos,dtype = np.float)
    print("Radius of the simulation:  "+str(rvir))
    # Logspace for radius
    logdelta = 3.0/50.0
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
            #m = 1
            #if i%m == 0:
                #title = lvl+"_"+halo+"_"+str(i//m)+".png"
                #plotHalo(title,a,b,c,posx,lvl,halo)
        else:
            axes = np.ones(3,dtype = np.float)*x_rad
            vecs = np.identity(3,dtype = np.float)

        rad = (abs(axes[0]*axes[1]*axes[2]))**(1./3.)
        #print("________________________________________________________________")
        print(x_rad,rad)
        if( x_rad > 2.0*rvir ):
            break
    semmiaxes.append(np.array([rvir,rvir,rvir]))
    semmiaxes = np.array(semmiaxes)
    np.savetxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt",semmiaxes, delimiter = ',')
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

