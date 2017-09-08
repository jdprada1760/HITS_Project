import matplotlib
matplotlib.use('Agg')

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

'''
#import yt
#from yt.units import parsec, Msun
from Py_Libs.shape import *

# C imports
import ctypes as ct
inertia = ct.cdll.LoadLibrary('./C_Libs/inertia.so')

# Simiulation specs
#lvl = 'level4_MHD'
lvl = 'level4_DM'
halonums = range(1,31)
snapnum = 127
#snapnum = 63
#snapnum = 255

# Keeps semiaxes and eigenvecs
semmiaxes = []
eigenvecs = []
for j in halonums:
    # Loads simulation
    print("--------------------------------------------------------------------------")
    print(j)
    halo = 'halo_'+str(j)
    pos,rvir = loadSim(lvl,halo,snapnum)
    pos = np.array(pos,dtype = np.float)
    print("Radius of the simulation:  "+str(rvir))
    # Pointers as parameters for C_function
    axes = np.ones(3,dtype = np.float)*rvir
    vecs = np.identity(3,dtype = np.float)
    # Calculates eigensystem of intertia tensor
    inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))
    semmiaxes.append(np.copy(axes))
    eigenvecs.append(np.copy(vecs))
    # Process is repeated for R = 0.01*rvir
    axes = np.ones(3,dtype = np.float)*0.01*rvir
    vecs = np.identity(3,dtype = np.float)
    inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))
    semmiaxes.append(np.copy(axes))
    eigenvecs.append(np.copy(vecs))

# Reshapes so that each halo info is in the same row
semmiaxes = np.reshape(np.array(semmiaxes),(len(halonums),6))
eigenvecs = np.reshape(np.array(eigenvecs),(len(halonums),18))
np.savetxt("../Plots/"+lvl+"/Semiaxes_rvir_1e2rvir.csv",semmiaxes,delimiter=',')
np.savetxt("../Plots/"+lvl+"/Eigenvecs_rvir_1e2rvir.csv",eigenvecs,delimiter=',')

''' 

# Simiulation specs
lvl = 'level4_MHD'
#lvl = 'level4_DM'
# Obtains the axes
axes = np.loadtxt("../Plots/"+lvl+"/Semiaxes_rvir_1e2rvir.csv", delimiter = ',')
# Axes in inner regions 0.01*rvir ~ 1-5 kpc
axes1 = axes[:,3:]
# Axes in outer regions rvir ~ 100-500 kpc
axes2 = axes[:,:3]

# Plots axial ratios c/a Vs b/a for R = 0.01Rvir 
plt.plot(axes1[:,1]/axes1[:,0],axes1[:,2]/axes1[:,0], marker = 'o',
 c = 'b', alpha = 0.6, linewidth = 0, label = r"$R = 0.01R_{vir}$" )
# Plots axial ratios c/a Vs b/a for R = Rvir 
plt.plot(axes2[:,1]/axes2[:,0],axes2[:,2]/axes2[:,0], marker = '^',
 c = 'r', alpha = 0.7, linewidth = 0, label = r"$R = R_{vir}$" )
# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'o')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0)#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '*',linewidth = 0)#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = '*',linewidth = 0)# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Triaxiality Inner-Outterskirts "+lvl)
plt.legend()
plt.savefig("../Plots/"+lvl+"/Triaxiality_"+lvl+".png")











