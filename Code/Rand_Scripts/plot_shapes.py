import numpy as np
import matplotlib
matplotlib.use('Agg')
from Py_Libs.shape import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

#lvl3 = 'level3_MHD'
#lvl4 = 'level4_MHD'
lvl3 = 'level3_DM'
lvl4 = 'level4_DM'
halonums = [6,16,21,23,24,27]
with PdfPages("../Plots/"+lvl3+"/Convergence.pdf") as pdf:
	for i in halonums:
	    halo = 'halo_'+str(i)
	    arr3 = np.loadtxt("../Plots/"+lvl3+"/"+halo+"/"+"abc_"+lvl3+"_"+halo+".txt", delimiter = ',')
	    arr4 = np.loadtxt("../Plots/"+lvl4+"/"+halo+"/"+"abc_"+lvl4+"_"+halo+".txt", delimiter = ',')
	    a3,b3,c3 = arr3[:-1].T
	    a4,b4,c4 = arr4[:-1].T
	    yvals3 = np.array([b3/a3,c3/a3,c3/b3])
	    xvals3 = (a3*b3*c3)**(1./3.)
	    yvals4 = np.array([b4/a4,c4/a4,c4/b4])
	    xvals4 = (a4*b4*c4)**(1./3.)
	    ylabel = ['b/a','c/a','c/b']
	    rvir = arr3[-1][0]
	    fig, axs = plt.subplots(figsize=(10,10),nrows=len(yvals3))
	    for ax,yval3,yval4,ylab in zip(axs,yvals3,yvals4,ylabel):
		 ax.plot(xvals3,yval3, c = 'b', label = "Level 3 MHD")
		 ax.plot(xvals4,yval4, c = 'r', label = "Level 4 MHD")
		 ax.plot([rvir,rvir],[0,1])
		 ax.set_xscale('log')
		 # Plotting ratios
		 ax.set_ylim(0,1)
		 # Valid for all Milkyway-like galaxies
		 ax.set_xlim(0.1,rvir+30)
		 ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		 ax.set_ylabel(ylab)
	    axs[-1].set_xlabel("log(R(kpc/h))")
	    pdf.savefig(fig)
	    plt.close()
