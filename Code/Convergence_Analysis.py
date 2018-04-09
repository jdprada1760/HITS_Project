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
with PdfPages("../Plots/Convergence_lvl34_"+sys.argv[1]+".pdf") as pdf:
    for j in halonums:
        # Minimum array at which there is enough particles for the shape to converge (among all graphics)
        mini = 0
        if sys.argv[1] == 'DM':
            mini = 121
        elif sys.argv[1] == 'MHD':
            mini = 148
        # Halo name        
        halo = "halo_"+str(j)        
    
        # Paths to files
        path3 = "../Plots/"+lvl3+"/"+halo+"/"+"abc_"+lvl3+"_"+halo+".txt"
        path4 = "../Plots/"+lvl4+"/"+halo+"/"+"abc_"+lvl4+"_"+halo+".txt"
    
        # Semiaxes
        arr3 = np.loadtxt(path3, delimiter = ",")[-mini:]
        arr4 = np.loadtxt(path4, delimiter = ",")[-mini:]
        rvir = arr3[-1][0]
        #rvir4 = arr4[-1][0]
        a3,b3,c3 = arr3[:-1].T
        a4,b4,c4 = arr4[:-1].T

        # Axial ratios and radii
        yvals3 = np.array([b3/a3,c3/a3,c3/b3])
        xvals3 = (a3*b3*c3)**(1./3.)
        yvals4 = np.array([b4/a4,c4/a4,c4/b4])
        xvals4 = (a4*b4*c4)**(1./3.)
         
        # The difference in axis ratios (we are omitting differences in a)
        diff = []
        count = 0.0

        # Finds the closest radius in level3 for each in level4
        for i4 in range(len(xvals4)):            
            
            # The index of closest radius in level3
            i3 = np.argmin(abs(xvals3-xvals4[i4]))

            # Appends the difference in axis ratios
            diff.append(abs((yvals3.T[i3]-yvals4.T[i4])/yvals3.T[i3]))

            if( abs(xvals3[i3] - xvals4[i4])/xvals3[i3] > 0.03 ):
                #print(halo+"____r= "+str(xvals3[i3])+"____delta= "+str(abs(xvals3[i3] - xvals4[i4])/xvals3[i3]))
                count += 1
        print(halo+"____Percent of points above 1% error in radius: " +str(count/len(xvals4)))
        diff = np.array(diff)
        
        # keeps track of standard deviation and mean
        #mini = 121
        xrand = np.zeros(mini)
        yrand = []
        
        # Fonts 
        MEDIUM_SIZE = 30
        SMALL_SIZE = 25
        SSSMALL_SIZE = 15
        plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

        # Graphics
        ylabel = ['b/a','c/a','c/b']
        fig, axs = plt.subplots(figsize=(10,10),nrows=len(yvals3))
        for ax,yval3,yval4,ylab,delta,indi in zip(axs,yvals3,yvals4,ylabel,diff.T,range(3)):
            
            yrand = []
            # Random sampling graphics
            for randn in range(10):
                path = "../Plots/"+lvl3+"/"+halo+"/rand_sample/"+str(randn)+".txt"
                arr = np.loadtxt(path, delimiter = ",")[-mini:]
                a,b,c = arr[:-1].T
                xvals = (a*b*c)**(1./3.)
                yvals = np.array([b/a,c/a,c/b])
                ax.plot(xvals,yvals[indi], c = 'grey', linewidth=0.1, zorder = 1)
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                # Statistics
                yrand.append(yvals[indi,-mini:])
                #print len(yvals[indi,-mini:])
                #print len(yvals[-mini:])
                if randn == 9:
                    xrand = xvals[-mini:]
                    
                # Hide axes
                plt.setp( ax.get_xticklabels(), visible=False)                
            
            yrand = np.array(yrand)
            
            #print(yrand)
            #print(xrand[-mini:])
            
            # Std and mean value
            #ax.plot(xrand[-mini:], np.mean(yrand,axis = 0)[-mini:], c = 'black', linewidth=0.2)
            
            yerr = np.zeros(len(xvals3))
            for l in range(len(yrand[0])):
                yerr[-(l+1)] = np.std(yrand[:,-(l+1)])
            
            #print yerr
            # Plots level3, level4 and their difference
            ax.errorbar(xvals3,yval3,yerr=0*yerr, c = 'blue', label = str(lvl3), linewidth=1.2, zorder = 2)
            
            ax.fill_between(xvals3, yval3 - 3*yerr,  
                 yval3 + 3*yerr, color="blue", alpha = 0.3, zorder = 1) 
                 
            ax.plot(xvals4,yval4, c = 'magenta', label = str(lvl4), linewidth=1.2, zorder = 10)
            ax.plot(xvals4,delta, c = 'g', label = "Difference", linewidth=1)
            ax.axhline(0.05, linewidth = 0.1)
            #ax.plot([rvir,rvir],[0,1], label = r"$R_{200}$")
            ax.plot([rvir,rvir],[0,1], label = r"$R_{200}$", linestyle = '--')
            ax.set_xscale('log')
            
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


            
            
            # Plotting ratios
            ax.set_ylim(0,1)
            
            # Axs specs
            ax.set_xlim(0.5,xvals3[-1]+30)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.set_ylabel(ylab)
            #ax.legend()
            
        # Axs specs    
        #plt.grid()
        axs[0].legend()
        axs[-1].set_xlabel("R(Kpc/h)")
        axs[-1].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axs[-1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.setp( axs[-1].get_xticklabels(), visible=True)                
        axs[0].set_title(halo)        
        pdf.savefig(fig)
        plt.close()



        
        
