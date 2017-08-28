import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Py_Libs.shape import *

# Simiulation specs
lvl = 'level4_MHD'
lvl2 = 'level4_DM'
halonums = range(1,3)

# Loads simulation
snapnum = 127

for j in halonums:
    print("--------------------------------------------------------------------------")
    print(j)
    halo = 'halo_'+str(j)
    pos1,rvir1 = loadSim(lvl,'halo_1',snapnum)
    pos1,rvir2 = loadSim(lvl2,'halo_1',snapnum)
    axes1 = np.loadtxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt", delimiter = ',')
    axes2 = np.loadtxt("../Plots/"+lvl2+"/"+halo+"/"+"abc_"+lvl2+"_"+halo+".txt", delimiter = ',')
    #axes11 = np.append(axes1,[[rvir1,rvir1,rvir1]], axis = 0)
    #axes22 = np.append(axes2,[[rvir2,rvir2,rvir2]], axis = 0)
    #np.savetxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt", axes11, delimiter = ',')
    #np.savetxt("../Plots/"+lvl2+"/"+halo+"/"+"abc_"+lvl2+"_"+halo+".txt", axes22, delimiter = ',')
    a,b,c = axes1[:-1].T
    rs = (a*b*c)**(1./3.)
    plt.plot(np.log10(rs),b/a)
    plt.plot(np.log10(rs),c/a)
    plt.savefig('tmp.png')

