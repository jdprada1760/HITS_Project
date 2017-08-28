import matplotlib
matplotlib.use('Agg')

from arepo import *
from shape import loadSim
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def fun1():
    lvlDM,lvlMHD = 'level4_DM','level4_MHD'
    lvl = 'level4_MHD'
    with PdfPages('../Plots/'+lvl+'/AxialRatios.pdf') as pdf:
        #arr = [1,2,3,4,5,7,8,10]
        arr = range(1,31)
        for j in arr:
            halo = 'halo_'+str(j)
            print(halo)
            path = "../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt"
            a,b,c = np.loadtxt(path,delimiter = ',').T
            snapnum = 127
            pos,rvir = loadSim(lvl,halo,snapnum)
            fig, ax1= plt.subplots(nrows=1, figsize=(15,7))
            yvals = np.array([b/a,c/a,c/b,(a**2-b**2)/(a**2-c**2)])
            ylabel = ['b/a','c/a','c/b','T']
            xvals = (a*b*c)**(1./3.)
            marks = ['-','-','-','--']
            #colors = ['r','g','b']
            for y,ylab,mark in zip(yvals,ylabel,marks):
                ax1.plot(xvals,y,mark,label =ylab )
            ax1.plot([rvir,rvir],[0,1])
            ax1.set_xscale('log')
            # Plotting ratios
            ax1.set_ylim(0,1.01)
            # Valid for all Milkyway-like galaxies
            ax1.set_xlim(0.1,2*rvir)
            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax1.set_xlabel("log(R(kpc/h))")
            ax1.legend()
            ax1.set_title(halo)
            pdf.savefig(fig)
            plt.close()


def fun2():
    lvlDM,lvlMHD = 'level4_DM','level4_MHD'
    print(lvlDM,lvlMHD)
    # lvl = 'level4_MHD'
    with PdfPages('../Plots/lvl4_DM_MHD_Comp.pdf') as pdf:
        #arr = [1,2,3,4,5,7,8,10]
        arr = range(1,31)
        for j in arr:
            halo = 'halo_'+str(j)
            print(halo)  
            snapnum = 127
            fig, axs= plt.subplots(nrows=4, figsize=(10,15))
            # Loads values
            pos,rvir1 = loadSim(lvlDM,halo,snapnum) 
            path = "../Plots/"+lvlDM+"/"+halo+"/"+"abc_"+lvlDM+"_"+halo+".txt"
            a,b,c = np.loadtxt(path,delimiter = ',').T
            yvals1 = np.array([b/a,c/a,c/b,(a**2-b**2)/(a**2-c**2)])
            xvals1 = (a*b*c)**(1./3.)
            # 
            pos,rvir2 = loadSim(lvlMHD,halo,snapnum)
            path = "../Plots/"+lvlMHD+"/"+halo+"/"+"abc_"+lvlMHD+"_"+halo+".txt"
            a,b,c = np.loadtxt(path,delimiter = ',').T
            yvals2 = np.array([b/a,c/a,c/b,(a**2-b**2)/(a**2-c**2)])
            xvals2 = (a*b*c)**(1./3.)
            # Graph Specs
            ylabel = ['b/a','c/a','c/b','T']
            #marks = ['-','-','-','--']
            #colors = ['r','g','b']
            for y1,y2,ylab,ax in zip(yvals1,yvals2,ylabel,axs):
                ax.plot(xvals1,y1,linestyle = '-',color ='g',label ="DM", linewidth = 2 )
                ax.plot(xvals2,y2,linestyle = '--',color = 'b',label ="MHD", linewidth = 1 )
                ax.plot([rvir1,rvir1],[0,1], linewidth = 2, linestyle = '-', color = 'g')
                ax.plot([rvir2,rvir2],[0,1], linewidth = 1, linestyle = '--', color = 'b')
                ax.set_ylabel(ylab)
                ax.set_xscale('log')
                # Plotting ratios
                ax.set_ylim(0,1.01)
                # Valid for all Milkyway-like galaxies
                ax.set_xlim(0.1,2*max([rvir1,rvir2]))
                ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                ax.legend()
            axs[-1].set_xlabel("log(R(kpc/h))")
            axs[0].set_title(halo)
            pdf.savefig(fig)
            plt.close()
