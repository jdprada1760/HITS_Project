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
        
        # Plots last snapshot in terms of radius
        arr = np.loadtxt(path+"/"+"abc_"+lvl+"_"+halo+".txt", delimiter = ",")
        a,b,c = (arr)[:-1].T

        # Plots    
        fig, axs = plt.subplots(figsize=(10,10))
        
        # Color settings
        my_norm = colors.Normalize(0,1)
        mapa = cm.ScalarMappable(norm=my_norm, cmap='hsv')
        
        # Loads redshift rvir shapes
        arr = np.loadtxt(path+"/"+"Z_axes.txt", delimiter = ",")
        redshift,a,b,c = arr.T        
        
        for j in range(len(a)):
            
            print("Redshift: "+str(redshift[j]))
            if redshift[j] > 1.5:
                continue
            
            # Plot scatter point
            my_col = mapa.to_rgba((1.0/(1.0+redshift[j])))
            axs.scatter(b[j]/a[j],c[j]/a[j], color = my_col, s=90, alpha = 1, zorder = len(a)-j)
            
            if redshift[j+1] < 1.5:
                axs.plot([b[j]/a[j],b[j+1]/a[j+1]],[c[j]/a[j],c[j+1]/a[j+1]], color = my_col, linewidth=2, alpha = 1, zorder = len(a)-j)
                
                
        # Plots last snapshot in terms of radius
        arr = np.loadtxt(path+"/"+"abc_"+lvl+"_"+halo+".txt", delimiter = ",")
        a,b,c = (arr)[:-1].T
        rvir =  arr[-1][0]
        xvals = (a*b*c)**(1./3.)
        # Axial ratios and radii (including rvir)
        a,b,c = (arr[xvals>=10])[:-1].T
        xvals = (a*b*c)**(1./3.)
        
        # Radius at which comparison is performed (in kpc)
        axs.plot((b/a)[xvals>rvir],(c/a)[xvals>rvir], color = 'black', linewidth=1.5, alpha = 1, zorder = len(a))
        axs.plot((b/a)[xvals<rvir],(c/a)[xvals<rvir], color = 'black', linewidth=1.5, alpha = 1, marker = "+", zorder = len(a))
    
        # Axs specs
        axs.plot([0,1],[0,1], linewidth= 2, c = 'black')
        axs.set_xlabel("b/a")
        axs.set_ylabel("c/a")
        axs.set_title("Triaxiality "+halo)
        axs.set_xlim(0,1)
        axs.set_ylim(0,1)
        #axs.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        # Colorbar
        xticks = [0, 0.5, 1, 2,4,8]
        yticks = [1.0/(1.0+x) for x in xticks]
        
        # Colorbar ax
        cbaxes = fig.add_axes([0.95, 0.1, 0.02, 0.8]) 
        cbar = mpl.colorbar.ColorbarBase(cbaxes, cmap = "hsv", norm = my_norm, orientation = 'vertical', ticks = yticks )  
        cbar.ax.set_yticklabels(xticks)
        cbar.set_label('$Redshift$', fontsize=14 )
        
        # Save
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
































