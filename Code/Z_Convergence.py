import matplotlib
matplotlib.use('Agg')

# Python
#from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import os
import sys

# Simiulation specs
lvl3 = sys.argv[1]
halonums = np.array([6,16,21,23,24,27])

#halonums = [24]
with PdfPages("../Plots/Redshift_"+lvl3+".pdf") as pdf:

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

        # Fonts 
        MEDIUM_SIZE = 30
        SMALL_SIZE = 25
        SSSMALL_SIZE = 17
        plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            
        # Plots    
        ylabel = ['b/a','c/a','c/b']
        fig, axs = plt.subplots(figsize=(10,15),nrows=3)
       
        # Color settings
        my_norm = colors.Normalize(0.0,0.8)
        mapa = cm.ScalarMappable(norm=my_norm, cmap="CMRmap")
        
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
                
                
                my_col = mapa.to_rgba(1.-(1.0/(1.0+redshift)))
                if(redshift < 1e-6):
                    mappable = ax.plot(xvals,yval, color = my_col, linewidth=1, label = "Radial shape", alpha = 1,zorder=10)
                    ax.plot([rvir,rvir],[0,0.5], color = my_col, linestyle = '-.',label = r"$R_{500}$",zorder =1)
                else:
                    mappable = ax.plot(xvals,yval, color = my_col, linewidth=0.5,  alpha = 1,zorder = 10)
                    ax.plot([rvir,rvir],[0,0.5], color = my_col,linestyle = '-.',zorder=1)
                
                # Major ticks every 20, minor ticks every 5
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
                ax.set_ylabel(ylab)
            
                # Plotting ratios
                ax.set_ylim(0,1)
        
                # Axs specs
                ax.set_xlim(3,rvir+30)
                ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                
                # Hide axes 
                plt.setp( ax.get_xticklabels(), visible=False)                
                
            
        # Axs specs  
        #axs[0].legend(loc='upper left', bbox_to_anchor=(1, 0.5))  
        
        xticks = [0, 0.5, 1, 2]
        yticks = [1.- 1.0/(1.0+x) for x in xticks]
        
        # Colorbar ax
        cbaxes = fig.add_axes([0.95, 0.1, 0.02, 0.8]) 
        cbar = mpl.colorbar.ColorbarBase(cbaxes, cmap = "CMRmap", norm = my_norm, orientation = 'vertical', ticks = yticks )  
        cbar.ax.set_yticklabels(xticks)
        cbar.set_label('$Redshift$', fontsize=30 )

        axs[-1].legend(loc = 0)
        axs[-1].set_xlabel("R(Kpc/h)")
        plt.setp( axs[-1].get_xticklabels(), visible=True)                
        #axs[0].set_title(halo+"  Rvir="+str(rvir)+"Kpc" )  
        axs[0].set_title(halo)  
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()


        
