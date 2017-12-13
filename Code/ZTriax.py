import matplotlib
matplotlib.use('Agg')

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import os

# Simiulation specs
lvl3 = 'level3_DM'
#lvl = 'level3_DM'
#halonums = range(1,31)
halonums = [6,16]

# Keeps semiaxes and eigenvecs
with PdfPages("../Plots/ZTriax_lvl3_DM.pdf") as pdf:

    for j in halonums:
    
        # Halo name        
        halo = "halo_"+str(j)     
        print("------------------------------------"+halo)   
    
        # Paths to files
        path = "../Plots/"+lvl3+"/"+halo
        
        # List of directories
        listdir = next(os.walk(path))[1]
        if 'rand_sample' in listdir:
            listdir.remove('rand_sample')
        def getkey(l):
            return int(l.split("_")[1])
        listdir.sort(key=getkey)

        # Plots    
        fig, axs = plt.subplots(figsize=(10,10))
        
        # Plots last snapshot in terms of radius
        arr = np.loadtxt(path+"/"+"abc_"+lvl3+"_"+halo+".txt", delimiter = ",")
        a,b,c = (arr)[:-1].T
        xvals = (a*b*c)**(1./3.)
        # Axial ratios and radii (including rvir)
        a,b,c = (arr[xvals>=9])[:-1].T
        
        # Radius at which comparison is performed (in kpc)
        radius = 50
        axs.plot(b/a,c/a, color = 'black', linewidth=1, alpha = 1)
       
        # Color settings
        my_norm = colors.Normalize(0,1)
        mapa = cm.ScalarMappable(norm=my_norm, cmap='hsv')
        
        for direc in listdir:
            
            # Semiaxes
            arr = np.loadtxt(path+"/"+direc+"/s"+direc.split("_")[1]+".txt", delimiter = ",")
            print(direc)
            
            # Loads information dictionary
            d = np.load(path+"/"+direc+"/info.npy")
            redshift = d.item().get('Redshift')
            print("Redshift: "+str(redshift))
            if redshift > 2:
                continue
     
            # Axial ratios and radii (including rvir)
            a,b,c = arr[:-1].T
            xvals = (a*b*c)**(1./3.)
            #yvals = np.array([b/a,c/a,c/b])
            rvir = arr[-1][0]
            
            # Finding the rradius which is nearer to the reference radius
            mini = np.argmin(abs(xvals-rvir))
            print("Best could do was : "+str(abs(xvals-radius)[mini]/radius))
            
            # Plot scatter point
            my_col = mapa.to_rgba((1.0/(1.0+redshift)))
            axs.scatter(b[mini]/a[mini],c[mini]/a[mini], color = my_col, s=90, alpha = 1)
    
        # Axs specs
        axs.plot([0,1],[0,1], linewidth= 2, c = 'black')
        axs.set_xlabel("b/a")
        axs.set_ylabel("c/a")
        axs.set_title("Triaxiality")
        axs.set_xlim(0,1)
        axs.set_ylim(0,1)
        #axs.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        # Colorbar
        xticks = [0, 0.5, 1, 2]
        yticks = [1.0/(1.0+x) for x in xticks]
        
        # Colorbar ax
        cbaxes = fig.add_axes([0.95, 0.1, 0.02, 0.8]) 
        cbar = mpl.colorbar.ColorbarBase(cbaxes, cmap = "hsv", norm = my_norm, orientation = 'vertical', ticks = yticks )  
        cbar.ax.set_yticklabels(xticks)
        cbar.set_label('$Redshift$', fontsize=14 )
        
        # Save
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
































