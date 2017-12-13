import matplotlib
matplotlib.use('Agg')

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import os


# Simiulation specs
lvl3 = "level3_DM"
halonums = np.array([6,16,21])

#halonums = [24]
with PdfPages("../Plots/Redshift_lvl3DM.pdf") as pdf:

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
        ylabel = ['b/a','c/a','c/b']
        fig, axs = plt.subplots(figsize=(10,10),nrows=3)
       
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
            yvals = np.array([b/a,c/a,c/b])
            rvir = arr[-1][0]
            mappable = 0
            for ax,yval,ylab in zip(axs,yvals,ylabel):
                
                
                my_col = mapa.to_rgba((1.0/(1.0+redshift)))
                mappable = ax.plot(xvals,yval, color = my_col, linewidth=0.5, label = "Redshift: "+ str(redshift), alpha = 1)
                ax.plot([rvir,rvir],[0,1], color = my_col)
                ax.set_xscale('log')
            
                # Plotting ratios
                ax.set_ylim(0,1)
        
                # Axs specs
                ax.set_xlim(3,rvir+20)
                ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                
            
        # Axs specs  
        #axs[0].legend(loc='upper left', bbox_to_anchor=(1, 0.5))  
        
        xticks = [0, 0.5, 1, 2]
        yticks = [1.0/(1.0+x) for x in xticks]
        
        # Colorbar ax
        cbaxes = fig.add_axes([0.95, 0.1, 0.02, 0.8]) 
        cbar = mpl.colorbar.ColorbarBase(cbaxes, cmap = "hsv", norm = my_norm, orientation = 'vertical', ticks = yticks )  
        cbar.ax.set_yticklabels(xticks)
        cbar.set_label('$Redshift$', fontsize=14 )

        
        axs[-1].set_xlabel("log(R(kpc/h))")
        axs[0].set_title(halo+"  Rvir="+str(rvir)+"Kpc" )    
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()


        
