import numpy as np
import matplotlib
matplotlib.use('Agg')
from Py_Libs.shape import *
import matplotlib.pyplot as plt

lvl = 'level4_MHD'
halonums = range(1,3)
for i in halonums:
    halo = 'halo_'+str(i)
    snapnum = 127  
    arr = np.loadtxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt", delimiter = ',')
    a,b,c = arr[:-1].T
    yvals = np.array([b/a,c/a,c/b])
    ylabel = ['b/a','c/a','c/b']
    xvals = (a*b*c)**(1./3.)
    rvir = arr[-1][0]
    fig,axs = ratiosPlots(xvals, yvals, ylabel, rvir)  
    plt.savefig("../Plots/"+lvl+"/"+halo+"/"+"tmp.png")
   
