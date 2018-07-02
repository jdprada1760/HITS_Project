import matplotlib
matplotlib.use('Agg')

# Python
#from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import os
import sys
# Simiulation specs
lvl = sys.argv[1]
#lvl = 'level3_DM'
#halonums = range(1,31)
halonums = [6,16,21,23,24,27]

# Keeps semiaxes and eigenvecs
with PdfPages("../Plots/"+lvl+"/"+"Z_Triax_"+lvl+".pdf") as pdf:

    for j in halonums:
    
        # Halo name        
        halo = "halo_"+str(j)     
        print("------------------------------------"+halo)   
    
        # Paths to file
        path = "../Plots/"+lvl+"/"+halo
        
        # Loads redshift rvir shapes
        arr = np.loadtxt(path+"/"+"Z_axes.txt", delimiter = ",")
        redshift,a,b,c = arr.T  

        # Fonts 
        MEDIUM_SIZE = 30
        SMALL_SIZE = 25
        SSSMALL_SIZE = 20
        plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        
        # Plots    
        fig, axs = plt.subplots(figsize=(10,10))
        
        # Color settings
        my_norm = colors.Normalize(0.0,0.8)
        mapa = cm.ScalarMappable(norm=my_norm, cmap='CMRmap')      
        
        for j in range(len(a)):
            
            
            if redshift[j] > 1.5:
                continue
            print("Redshift: "+str(redshift[j]))
            # Plot scatter point
            my_col = mapa.to_rgba(1-(1.0/(1.0+redshift[j])))
            if redshift[j] < 1e-6:
            
                axs.scatter(b[j]/a[j],c[j]/a[j], color = my_col, s=90, alpha = 1, zorder = len(a)-j,label = r"Historical shape $R_{500}$")
            else:
                axs.scatter(b[j]/a[j],c[j]/a[j], color = my_col, s=90, alpha = 1, zorder = len(a)-j)
            
            if redshift[j+1] < 1.5:
                axs.plot([b[j]/a[j],b[j+1]/a[j+1]],[c[j]/a[j],c[j+1]/a[j+1]], color = my_col, linewidth=3, alpha = 1, zorder = len(a)-j)
            


                
                
        # Plots last snapshot in terms of radius
        arr = np.loadtxt(path+"/"+"abc_"+lvl+"_"+halo+".txt", delimiter = ",")
        a,b,c = (arr)[:-1].T
        rvir =  arr[-1][0]
        xvals = (a*b*c)**(1./3.)
        # Axial ratios and radii (including rvir)
        a,b,c = (arr[:-1][xvals>=1]).T
        xvals = (a*b*c)**(1./3.)
        
        # Radius at which comparison is performed (in kpc)
        axs.plot((b/a)[xvals>rvir],(c/a)[xvals>rvir], color = 'black', linewidth=2, alpha = 1, linestyle = "--", zorder = len(a))
        axs.plot((b/a)[xvals<rvir],(c/a)[xvals<rvir], color = 'black', linewidth=2, alpha = 1, linestyle = "--", zorder = len(a),label="Radial shape z = 0")
    
        # Axs specs
        axs.plot([0,1],[0,1], linewidth= 2, c = 'black')
        axs.set_xlabel("b/a")
        axs.set_ylabel("c/a")
        axs.set_title(halo)
        
         # grid!!
        # Major ticks every 20, minor ticks every 5
        major_ticksy = np.linspace(0, 1, 3)
        minor_ticksy = np.linspace(0, 1, 11)
        #major_ticksx = np.logspace(-1, 2, 3)
        #minor_ticksx = np.logspace(-1, 2, 30)

        axs.set_xticks(major_ticksy)
        axs.set_xticks(minor_ticksy, minor=True)
        axs.set_yticks(major_ticksy)
        axs.set_yticks(minor_ticksy, minor=True)
        axs.grid(which='both')
        axs.grid(which='minor', alpha=0.3)
        axs.grid(which='major', alpha=0.6)

        axs.set_xlim(0,1)
        axs.set_ylim(0,1)
        #axs.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        # Colorbar
        xticks = [0, 0.5, 1,2]
        yticks = [1.- 1.0/(1.0+x) for x in xticks]
        
        # Colorbar ax
        cbaxes = fig.add_axes([0.95, 0.1, 0.02, 0.8]) 
        
        axs.legend(loc=0)
        cbar = mpl.colorbar.ColorbarBase(cbaxes, cmap = "CMRmap", norm = my_norm, orientation = 'vertical', ticks = yticks )  
        cbar.ax.set_yticklabels(xticks)
        cbar.set_label('$Redshift$', fontsize=30 )
        
        # Save
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
































