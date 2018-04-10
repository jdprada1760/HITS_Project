import matplotlib
matplotlib.use('Agg')
import sys

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages

# Simiulation specs
lvl3 = "level3_" +sys.argv[1]
lvl4 = "level4_" +sys.argv[1]
halonums = np.array([6,16,21,23,24,27])
#halonums = [24]
with PdfPages("../Plots/Simple_Convergence_lvl34_"+sys.argv[1]+".pdf") as pdf:
    for j in halonums:
        
        # Halo name        
        halo = "halo_"+str(j)        
    
        # Paths to files
        path3 = "../Plots/"+lvl3+"/"+halo+"/"+"abc_"+lvl3+"_"+halo+".txt"
        path4 = "../Plots/"+lvl4+"/"+halo+"/"+"abc_"+lvl4+"_"+halo+".txt"
    
        # Semiaxes
        arr3 = np.loadtxt(path3, delimiter = ",")
        arr4 = np.loadtxt(path4, delimiter = ",")
        rvir = arr3[-1][0]
        #rvir4 = arr4[-1][0]
        a3,b3,c3 = arr3[:-1].T
        a4,b4,c4 = arr4[:-1].T

        # Axial ratios and radii
        yvals3 = np.array([b3/a3,c3/a3,c3/b3])
        xvals3 = (a3*b3*c3)**(1./3.)
        yvals4 = np.array([b4/a4,c4/a4,c4/b4])
        xvals4 = (a4*b4*c4)**(1./3.)
        
        # Fonts 
        MEDIUM_SIZE = 30
        SMALL_SIZE = 25
        SSSMALL_SIZE = 17
        plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

        # Graphics
        ylabel = ['b/a','c/a','c/b']
        fig, axs = plt.subplots(figsize=(10,10),nrows=len(yvals3))
        for ax,yval3,yval4,ylab,indi in zip(axs,yvals3,yvals4,ylabel,range(3)):
            
            # Hide axes
            plt.setp( ax.get_xticklabels(), visible=False)                
                        
            #print(yrand)
            #print(xrand[-mini:])
            
            # Std and mean value
            #ax.plot(xrand[-mini:], np.mean(yrand,axis = 0)[-mini:], c = 'black', linewidth=0.2)
            
            
            #print yerr
            # Plots level3, level4 and their difference
                 
            ax.plot(xvals4,yval4, c = 'magenta', label = "level4_"+sys.argv[1], linewidth=2)
            ax.plot(xvals3,yval3, c = 'g', label = "level3_"+sys.argv[1], linewidth=2)
            #ax.axhline(0.05, linewidth = 0.1)
            #ax.plot([rvir,rvir],[0,1], label = r"$R_{200}$")
            ax.plot([rvir,rvir],[0,1], label = r"$R_{200}$", linestyle = '--')
            
            
            # Major ticks every 20, minor ticks every 5
            # grid!!
            major_ticksy = np.linspace(0, 1, 3)
            minor_ticksy = np.linspace(0, 1, 11)
            major_ticksx = np.logspace(-1, 2, 3)
            minor_ticksx = np.logspace(-1, 2, 30)

            ax.set_xticks(major_ticksx)
            ax.set_xticks(minor_ticksx, minor=True)
            ax.set_yticks(major_ticksy)
            ax.set_yticks(minor_ticksy, minor=True)
            ax.grid(which='both')
            ax.grid(which='minor', alpha=0.3)
            ax.grid(which='major', alpha=0.6)

            ax.set_xscale('log')
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            
            # Plotting ratios
            ax.set_ylim(0,1)
            
            # Axs specs
            ax.set_xlim(0.1,xvals3[-1]+30)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.set_ylabel(ylab)
            #ax.legend()
            
        # Axs specs    
        #plt.grid()
        axs[-1].legend(loc = 0)
        axs[-1].set_xlabel("R(Kpc/h)")
        axs[-1].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axs[-1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.setp( axs[-1].get_xticklabels(), visible=True)                
        axs[0].set_title(halo)        
        pdf.savefig(fig)
        plt.close()
