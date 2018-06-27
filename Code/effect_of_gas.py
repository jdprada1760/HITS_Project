########################################################################################################
##### THIS SCRITP IS FOR READING THE SUMARY OF IMPORTANT QUANTITIES FOR DM HALO SHAPE               ####
##### READS:                                                                                        ####
#####               PARAMETERS .CSV (GAS,STELLAR,BH MASSES, SCALE RADII & OTHER IMPORTANT PARAMS)   ####
#####               AXES & VECS.CSV (TRIAXIAL CHARACTERIZATION AT 4to5 DIFFERENT RADII)             ####
#####                               MEAN (VOLUMETRIC) DENSITY CHARACTERIZATION                      ####
#####                               ISOPOTENTIAL SHELL CHARACTERIZATION                             ####
#####                               ISODENSITY SHELL CHARACTERIZATION                               ####
#####                                                                                               ####
##### FILENAMES ARE INSIDE THE SCRIPT                                                               ####
########################################################################################################

import numpy as np
import csv
import os

#########################################################################################################

lvl = 'level4_'# level of the simullation
nhalos = 30 # number of halos

#lvl = 'level3_'# level of the simullation
#nhalos = 6 # number of halos

folder = '../Plots/' # The folder containing all files

# FILENAMES
# Contains parameters of the disk (MHD)
fn_params = lvl+'params'+'.csv'
# Triaxial properties different radii (MHD & DM) using volumetric (mean) density: Allgood et al.
fn_axes = lvl+'axes'+'.csv' # Axes = axial ratios (sqrt of eigenvals)
fn_vecs = lvl+'vecs'+'.csv' # Vecs = principal axes directions
# Triaxial properties using a simple approach to measure the isopotential shape
fn_potaxes = lvl+'potential_axes'+'.csv'
fn_potvecs = lvl+'potential_vecs'+'.csv'
# Triaxial properties using a simple approach to measure isodensity shells
fn_denaxes = lvl+'density_axes'+'.csv'
fn_denvecs = lvl+'density_vecs'+'.csv'


########################################################################################################

def main():

    ##############
    # Parameters #
    ##############
    params_nms,params = read_csv(folder+fn_params)
    # skip the first row of parameters (haloname)
    params = params.T[1:].T
    params_nms = params_nms[1:]
    
    # Dictionary of parameters
    par = dict(zip(params_nms,params.T))
    
    # Column fields are:
    #'SFR', 'DiskRad', 'MassGas', 'MassDM', 'MassStars', 'MassBH', 'HaloRadCrit500', 'vDisk1', 'vDisk2', 'vDisk3'
    #print(par,params_nms)
    
    # Important quantities (normalized)
    GasDiskRad = par['GasDiskRad']               # halfmass radius of gas
    StarDiskRad = par['StellarDiskRad']          # Same for stellar mass     
    SFRD = par['SFR']/(StarDiskRad**2)           # Star formation ratio density
    gas_den = par['MassGas']/(GasDiskRad**2)     # Mass densities
    star_den = par['MassStars']/(StarDiskRad**2)
    BH_den = par['MassBH']/(StarDiskRad**2)
    bar_frac = (par['MassGas']+par['MassStars'])/par['MassDM'] # The baryonic fraction
    
    # The disk directions 
    disk_vec = params.T[-3:].T

    ################
    # Axial_Ratios #
    ################
    tmp,axes = read_csv(folder+fn_axes)
    
    #################################################################################
    # Axes are sampled at 0.125R, 0.25R, 0.5R, R, RDisk         R is Rvir(Rcrit500)
    # axes[halo_index, 0, : ] is the vector (a,b,c) at 0.125R 
    #################################################################################
    
    axes = np.reshape(axes,(2*nhalos,5,3)) # contains DM and MHD (2*nhalos) axes (3) sampled at (5) radii
    axesMHD = axes[::2] #organized halo1MHD,halo1DM,halo2MHD,halo2DM....
    axesDM = axes[1::2]
    
    # Isopotential and Isodensity axes
    tmp,axes = read_csv(folder+fn_denaxes, skip_header = 0) # Isodensity (den)
    axesDen = np.reshape(axes,(nhalos,5,3)) # contains only MHD (nhalos) axes (3) sampled at (5) radii
    tmp,axes = read_csv(folder+fn_potaxes, skip_header = 0) # Isopotential (pot)
    axesPot = np.reshape(axes,(nhalos,5,3))
    
    #print axesMHD
    #print axesDM
     
    # Define the triaxiality parameter T 
    T_MHD = get_T(axesMHD) #shape: (nhalos,5)
    T_DM  = get_T(axesDM)
    
    # Defines Asphericity
    asph_MHD = asphericity(axesMHD) #shape: (nhalos,5)
    asph_DM  = asphericity(axesDM)
    
    # Define deviation from sphericity 
    
    
    ####################
    # Axial_Directions #
    ####################
    tmp,vecs = read_csv(folder+fn_vecs)
    
    # Vecs are sampled at 0.125R, 0.25R, 0.5R, R, RDisk         R is Rvir(Rcrit500)
    #################################################################################
    # Vecs are sampled at 0.125R, 0.25R, 0.5R, R, RDisk         R is Rvir(Rcrit500)
    # vecs[halo_index, 0, : ] is the matrix [va,vb,vc]  at 0.125R 
    #################################################################################
    
    vecs = np.reshape(vecs,(2*nhalos,5,3,3)) # contains DM and MHD (2*nhalos) vecs (3,3) sampled at (5) radii
    vecsMHD = vecs[::2]
    vecsDM = vecs[1::2]
    
    # Isopotential and Isodensity axes
    tmp,vecs = read_csv(folder+fn_denvecs, skip_header = 0) #skipHeader:Little error of unified notation (no problem,to fix later)
    vecsDen = np.reshape(vecs,(nhalos,5,3,3))
    tmp,vecs = read_csv(folder+fn_potvecs, skip_header = 0)
    vecsPot = np.reshape(vecs,(nhalos,5,3,3))

    ####################
    # Some graphics    #
    ####################
    quant = T_MHD - T_DM# Some measure of symmetry
    # np.abs(2*(0.5 - T_DM)) is 1 when object is axisymetric and 0 when it is maximally triaxial
    name = 'Delta Triaxiality'
    
    
    import matplotlib.pyplot as plt
    plt.plot(axesMHD[:,0,1]/axesMHD[:,0,0],axesMHD[:,0,2]/axesMHD[:,0,0], 'ro' )
    plt.plot(axesDM[:,0,1]/axesDM[:,0,0],axesDM[:,0,2]/axesDM[:,0,0], 'bo' )
    plt.show()
    
    
    simple_plot(GasDiskRad,quant[:,0],'GasDiskRad',name, name+' Vs GasDiskRad')
    simple_plot(StarDiskRad,quant[:,0],'StarDiskRad',name, name+' Vs StarDiskRad')
    simple_plot(SFRD,quant[:,0],'SFR density',name, name+' Vs SFR density', logx =True)
    simple_plot(gas_den,quant[:,0],'Gas Density',name, name+' Vs Gas Density', logx =True)
    simple_plot(star_den,quant[:,0],'Star Density',name, name +' Vs Star Density', logx =True)
    simple_plot(BH_den,quant[:,0],'BH Density',name, name+' Vs BH Density', logx =True)
    simple_plot(bar_frac,quant[:,0],'Baryonic fraction',name, name+' Vs Baryonic fraction')  
    
    
    quant = asph_MHD - asph_DM
    name = 'Delta asphericity'
    simple_plot(GasDiskRad,quant[:,0],'GasDiskRad',name, name+' Vs GasDiskRad')
    simple_plot(StarDiskRad,quant[:,0],'StarDiskRad',name, name+' Vs StarDiskRad')
    simple_plot(SFRD,quant[:,0],'SFR density',name, name+' Vs SFR density', logx =True)
    simple_plot(gas_den,quant[:,0],'Gas Density',name, name+' Vs Gas Density', logx =True)
    simple_plot(star_den,quant[:,0],'Star Density',name, name +' Vs Star Density', logx =True)
    simple_plot(BH_den,quant[:,0],'BH Density',name, name+' Vs BH Density', logx =True)
    simple_plot(bar_frac,quant[:,0],'Baryonic fraction',name, name+' Vs Baryonic fraction') 
    
    
    # We Want to study the alignment between the different vectorss
    radii = ['12','25','50','1','Disk']
    radii = ['Rad_'+ra for ra in radii]
    
    from matplotlib.backends.backend_pdf import PdfPages
    
    with PdfPages("./pics/alignments"+".pdf") as pdf:
        for j in range(30):
            # Graphics
            polar_graphs(pdf,vecsMHD[j],disk_vec[j],radii,"halo_"+str(j))
        
        
            # Text 
            print "______________________________________________________________"
            print "____________" + "halo_"+str(j)+"______________________________"
            print "______________________________________________________________\n"
            
            for i in range(5):
                print "############        " + radii[i] +"             ###############\n\n"
                print "b/a = " + str(axesMHD[j,i,1]/axesMHD[j,i,0]) 
                print "c/a = " + str(axesMHD[j,i,2]/axesMHD[j,i,0]) + "\n"
                
                print "Angle with Major axis"
                print (180.0/np.pi)*np.arccos(np.abs(np.sum(vecsMHD[j,i,0]*disk_vec[j]))) # both are normalized vectors
                #print vecsMHD[j,i,0], np.linalg.norm(vecsMHD[j,i,0])
                #print disk_vec[j], np.linalg.norm(disk_vec[j])
                print "\n"
                print "Mid axis"
                print (180.0/np.pi)*np.arccos(np.abs(np.sum(vecsMHD[j,i,1]*disk_vec[j]))) # both are normalized vectors
                print "\n"
                print "Minor axis"
                print (180.0/np.pi)*np.arccos(np.abs(np.sum(vecsMHD[j,i,2]*disk_vec[j]))) # both are normalized vectors
                print "\n"


def simple_plot(x,y,xlabel,ylabel,title, mode = 'scatter', logx = False):
    
    import matplotlib.pyplot as plt
    
    marker = '-'
    linestyle = '-'
    if mode == 'scatter' :
        marker = 'o'
        linestyle = 'None'
    if logx == False:
        plt.plot(x,y, marker = marker, linestyle = linestyle)
    else:
        plt.semilogx(x,y, marker = marker, linestyle = linestyle)
        
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    #plt.show()
    
    # Tries to create directory
    path = './pics/'
    try:
        os.makedirs(path)
    except OSError:
        print("Directory already exists")
    plt.savefig(path+title)
    plt.close()

# Reads file in "path/filename"
# Returns data as a numpy array
# Returns filenames as a list of strings

def read_csv(filename, skip_header = 1):

    # Removes quotations from file (little writing error) 
    os.system('sed -i \'s/"//g\' '+filename) 
    
    names = []
    with open( filename , 'r') as myfile:
        rd = csv.reader(myfile)
        names = next(rd) # gets only first line 
    
    data = np.genfromtxt(filename, delimiter = ',', skip_header =skip_header, dtype=float )
    
    return names,data

# Triaxialities T = (1-(b/a)^2)/(1-(c/a)^2) 
# (T = 0 --> Oblate ) (T = 1 --> Prolate ) (T = 0.5 --> Max Triax )
def get_T(axes):
   return (axes[:,:,0]**2-axes[:,:,1]**2)/(axes[:,:,0]**2-axes[:,:,2]**2)               
        
# Calculates a measure of asphericity
# Returns the 1 minus actual distance from the point in the plane (c/a,b/a) to the point (1,1) in the triaxiality plane
# asphericity is 1 if object is a sphere and bigger than 0 if axial ratios are different from 1
def asphericity(axes):
    return 1-(1./np.sqrt(2))*np.sqrt((1-(axes[:,:,1]/axes[:,:,0]))**2 + (1-(axes[:,:,2]/axes[:,:,0]))**2)    
        
# Given a vector matrix mat, and a sigle vector vec, rotates vector and returns its representation in coordinates theta,phi
# theta,phi = 0,0 : minor axis (z)
# theta,phi = 90,0: medium axis (x)
# theta,phi = 90,90: major axis (y)
# returns theta,phi (clearly)

def sphere_coords(mat,vec):
    
    # we asume mat = [vmajor,vmedium,vminor]
    z = np.abs(np.sum(mat[2]*vec)) # abs to avoid angles bigger than 90 deg
    y = np.sum(mat[0]*vec)
    x = np.sum(mat[1]*vec)
    
    # To spherical coordinates
    theta = (180.0/np.pi)*np.arccos(z)
    phi = (180.0/np.pi)*np.arctan(y/x)
    
    return theta,phi
    
# Given pdf makes, for a single halo makes a graphic of alignment of axes    
# Parameters: pdf              : The pdf manager to save fig
#             vecs(5,3,3)      : The vectors of reference at each radius  
#             disk_vec(3),     : The disk vector, unchanged in principle, but with different representations (calculated within)
#             radii            : The radii at which each triaxiality is measured      
def polar_graphs(pdf,vecs,disk_vec,radii,title):        
    
    import matplotlib.pyplot as plt
    
    # to modify the size
    fig = plt.figure(figsize=(15,10))
    fig.suptitle(title, fontsize=16)
    #  theta and phi
    ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2, projection = 'polar')
    ax2 = plt.subplot2grid((2,6), (0,2), colspan=2, projection = 'polar')
    ax3 = plt.subplot2grid((2,6), (0,4), colspan=2, projection = 'polar')
    ax4 = plt.subplot2grid((2,6), (1,1), colspan=2, projection = 'polar')
    ax5 = plt.subplot2grid((2,6), (1,3), colspan=2, projection = 'polar')
    
    axs = [ax1,ax2,ax3,ax4,ax5]
    for i in range(5):
        
        theta,phi = sphere_coords(vecs[i],disk_vec ) # Gets disk in the spherical coordinates determined by principal axes
        axs[i].set_rmax(110)
        axs[i].plot((np.pi/180.)*phi,theta, marker = 'o', c = 'black') # plot disk
        axs[i].set_rticks([30,45,60,90,110])
        axs[i].grid(True)
        axs[i].plot(90,0, marker = '*', c = 'y', markersize = 20) # plot minor axis
        axs[i].plot(0,90, marker = 's', c = 'b', markersize = 20) # plot mid axis
        axs[i].plot((np.pi/180.)*90,90, marker = '^', c = 'g', markersize = 20) # plot major axis
        #axs[i].plot((np.pi/180.)*110,0, linewidth = 0)#
        axs[i].set_title(radii[i])
        axs[i].margins(0)
        #axs[i].set_rlabel_position(95)
        #axs[i].legend()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()    
        
        



main()





















