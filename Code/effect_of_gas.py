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
    tmp,axes = read_csv(folder+fn_denaxes) # Isodensity (den)
    axesDen = np.reshape(axes,(nhalos,5,3)) # contains only MHD (nhalos) axes (3) sampled at (5) radii
    tmp,axes = read_csv(folder+fn_potaxes) # Isopotential (pot)
    axesPot = np.reshape(axes,(nhalos,5,3))
    
    #print axesMHD
    #print axesDM
     
    # Define the triaxiality parameter T 
    T_MHD = get_T(axesMHD) #shape: (nhalos,5)
    T_DM  = get_T(axesDM)
    
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
    tmp,vecs = read_csv(folder+fn_denvecs) 
    vecsDen = np.reshape(vecs,(nhalos,5,3,3))
    tmp,vecs = read_csv(folder+fn_potvecs)
    vecsPot = np.reshape(vecs,(nhalos,5,3,3))

    ####################
    # Lasso Fit        #
    ####################
    
    
    # Defines Sphericity
    asph_MHD = asphericity(axesMHD) #shape: (nhalos,5)
    asph_DM  = asphericity(axesDM)
    y = (asph_MHD - asph_DM)[:,0] # 0 -> only calculated at 0.125 Rvir
    
    #y1 = axesMHD[:,0,1]/axesMHD[:,0,0] # this is b/a MHD for all haloes
    #y2 = axesMHD[:,0,2]/axesMHD[:,0,0] # c/a MHD
    
    # The X values (Parameters of the fit)
    
    # Just to see if the increase is related to the unperturbed state of the halo
    #x1 = axesDM[:,0,1]/axesDM[:,0,0] #b/a DM (q)
    #x2 = axesDM[:,0,2]/axesDM[:,0,0] #c/a MHD (s)
    X = np.array([GasDiskRad,StarDiskRad,np.log10(SFRD),np.log10(gas_den),np.log10(star_den),np.log10(BH_den),bar_frac]).T
        
    Xnames = ['Gas Radius','StarDisk Radius','SFR density','Gas density','Star density','BH density',
                'Baryonic fraction']
    Xlog = [0,0,1,1,1,1,0,0,0,0] # To know if scales are logarithmic


    
    LASSO_fit(X,y,Xnames=Xnames, yname = 'Delta sphericity', Xlog = Xlog)
    #LASSO_fit(np.array([x1]).T,y,['q'])
    


    

# REturns the LASSO fit according the given yparams and xparams
# Plots important projections over the x space to see accuracy of fit.
# Returns parameters of the fit.

def LASSO_fit(X,y,Xnames='',yname='',Xlog = None):        
    
    import matplotlib.pyplot as plt
    # Import scikit learn
    from sklearn import linear_model
    # We apply lasso with cross validation to choose an optimal weight for penalty (alpha) (Not good results)
    clf = linear_model.LassoCV(normalize = True)
    
    clf = linear_model.Lasso(alpha = 0.002, normalize = True)
    
    # Asuming array x is of size (npoints,mparams) and array y is of size (npoints)
    clf.fit(X,y)
    
    Xs = np.array([np.linspace(min(x)-0.01,max(x)+0.01,100) for x in X.T])
    #print Xs
    
    ypredict = clf.predict(X)
    for i in range(len(Xnames)):
        print Xnames[i],
        print '\n',
        print clf.coef_[i],'\n'
        
        '''
        plt.plot(X[:,i],y,'bo')
        xmin = min(X[:,i])-0.05
        xmax = max(X[:,i])+0.05
        Xd = np.linspace(xmin,xmax,100)
        plt.plot(X[:,i],ypredict,'ro')
        plt.xlabel(Xnames[i])
        plt.show()
        plt.close() 
        '''
        
      
    print clf.intercept_      
    
    
    
    
    
    
def simple_plot(x,y,xlabel,ylabel,title, mode = 'scatter', logx = False, folder = ''):
    
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
    path = './pics/'+lvl.split('_')[0]+'/'
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
        
# Given a vector matrix mat, and a sigle vector vec, rotates vector and returns its representation (in mat) in coordinates theta,phi
# theta,phi = 0,0 : minor axis (z)
# theta,phi = 90,0: medium axis (x)
# theta,phi = 90,90: major axis (y)
# ref      = axis of reference to solve sign degeneracy
# reference vector is represented in the basis of Mat
# this reference should be  the previous value of the vector
# returns theta,phi (clearly)

def sphere_coords(mat,vec, ref = None):
    
    # we asume mat = [vmajor,vmedium,vminor]
    z = np.dot(mat[2],vec) # (abs) to avoid angles bigger than 90 deg
    y = np.dot(mat[0],vec)
    x = np.dot(mat[1],vec)
    
    posi = np.array([x,y,z])# Position vector

    
    # If no reference is given then we take the z axis as vector of reference
    if ref is None:
        if posi[2] < 0 :
            posi = -posi
            
    
    # If there is a vector of reference, then we take the sign which gives the smallest angle
    else:
        
        # First we check that there is no flip, the threshold is arbitrary    
        if np.arccos(np.abs(np.dot(posi,ref))) > 45:
            ref = np.ones(3)/np.sqrt(3) #If ther is a flipp, choose the mean vector of axes.
    
        #print "posi          ",posi
        sign = np.dot(posi,ref)#The sign along the reference vector
        #print "dot           ",sign
        #print ref
        if sign < 0 :
            posi = -posi
        #print "posi after    ",posi
        
    
    # To spherical coordinates
    theta = (180.0/np.pi)*np.arccos(posi[2])
    phi = (180.0/np.pi)*np.arctan(posi[1]/posi[0])
    if posi[0] < 0:
        phi += 180
    #print "theta,phi     ",theta,phi
    #print '--------'
    return theta,phi

# Returns the corresponding vector in cartesians
# theta is in degrees, phi is in radians (matplotlib )
def to_Cartesians(theta,phi):
    theta = (np.pi/180.)*theta
    return np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta) ])
    
    
    
# Given pdf makes, for a single halo makes a graphic of alignment of axes    
# Parameters: pdf              : The pdf manager to save fig
#             vecs(5,3,3)      : The vectors of reference at each radius  
#             disk_vec(3),     : The disk vector, unchanged in principle, but with different representations (calculated within)
#             radii            : The radii at which each triaxiality is measured      
def polar_graphs(pdf,vecs,disk_vec,radii,title, radii_vals = None, mode = 'all_in_1'):        
    
    import matplotlib.pyplot as plt
    
    fig =  None
    theta,phi = sphere_coords(vecs[0],disk_vec ) # Gets disk in the spherical coordinates determined by principal axes
    phi = (np.pi/180.)*phi # phi input must be in radians
    
    ######################################################################################################################
    ######################################################################################################################
    if mode == 'all_in_1':
        # to modify the size
        fig = plt.figure(figsize=(8,8))

    
        ax = fig.add_subplot(111, projection='polar')
        ax.set_rmin(0)
        ax.set_rmax(100)
        ax.set_rlim(0,100)
        ax.set_rlabel_position(-22.5)  
        ax.set_rticks([30,45,60,90,100])
        ax.set_axisbelow(True)
        ax.grid(True) 
        #ax.set_title(radii[i])
        ax.scatter([phi],[theta], marker = 'X', c = 'black', s = 250) # plot disk
        
        mint,minp = 0,90 # angles theta,phi of first radius (12,4%Rvir)
        midt,midp = 90,0
        majt,majp = 90,90
        minp = (np.pi/180.)*minp # phi must be in radians
        midp = (np.pi/180.)*midp
        majp = (np.pi/180.)*majp
        
        # The colormap # Color settings
        import matplotlib.cm as cm
        import matplotlib.colors as colors
        
        my_norm = colors.Normalize(np.log(0.8*min(radii_vals)),np.log(1.2*max(radii_vals)))
        mapa = cm.ScalarMappable(norm=my_norm, cmap='inferno')
        
        my_col = mapa.to_rgba(np.log(radii_vals[0]))
        # Axes of reference
        ec = ['green','red']# EdgeColors         
        ax.scatter([90],[0], marker = '*', c = my_col, s = 450, alpha = 0.7,
        linewidths = 1.5, edgecolors = 'green' ) # plot minor axis
        ax.scatter([0],[90], marker = 's', c = my_col, s = 300, alpha = 0.7,
        linewidths = 1.5, edgecolors = 'green') # plot mid axis
        ax.scatter([(np.pi/2)],[90], marker = '^', c = my_col, s =300, alpha = 0.7,
        linewidths = 1.5, edgecolors = 'green') # plot major axis

        for i in range(1,5):
        
            my_col = mapa.to_rgba(np.log(radii_vals[i]))
        
            mint,minp = sphere_coords(vecs[0,:,:],vecs[i,2,:],to_Cartesians(mint,minp))# angular parameters min
            midt,midp = sphere_coords(vecs[0,:,:],vecs[i,1,:],to_Cartesians(midt,midp))# angular parameters mid
            majt,majp = sphere_coords(vecs[0,:,:],vecs[i,0,:],to_Cartesians(majt,majp))# angular parameters maj
            
            minp = (np.pi/180.)*minp# radians for matplotlib
            midp = (np.pi/180.)*midp
            majp = (np.pi/180.)*majp
            
            '''
            print "Minor: "
            print vecs[0,[1,0,2],:].dot(vecs[i,2,:])
            print to_Cartesians(mint,minp)
            print "Mid: "
            print vecs[0,[1,0,2],:].dot(vecs[i,1,:])
            print to_Cartesians(midt,midp)
            print "Major: "
            print vecs[0,[1,0,2],:].dot(vecs[i,0,:])
            print to_Cartesians(majt,majp)
            print "____________________________________________________"
            '''
            
            #print mint, ec[int(mint>90)], int(mint>90), int(mint<90), ec
            ax.scatter([minp],[90-np.abs(90-mint)], marker = '*', c = my_col,s = 450, alpha = 0.7, 
            linewidths = 1.5, edgecolors = ec[int(mint>90)] )
            ax.scatter([midp],[90-np.abs(90-midt)], marker = 's', c = my_col,s = 300, alpha = 0.7, 
            linewidths = 1.5, edgecolors = ec[int(midt>90)])
            ax.scatter([majp],[90-np.abs(90-majt)], marker = '^', c = my_col,s = 300, alpha = 0.7, 
            linewidths = 1.5, edgecolors = ec[int(majt>90)])
            
        # Colorbar ax
        import matplotlib as mpl
        from matplotlib.ticker import FormatStrFormatter
        cbaxes = fig.add_axes([0.95, 0.1, 0.02, 0.8])
        
        xticks = np.array([0.10, 0.20, 0.40, 0.80, 1])*radii_vals[-2]
        yticks = np.log(xticks)  
        cbar = mpl.colorbar.ColorbarBase(cbaxes, cmap = "inferno", norm = my_norm, orientation = 'vertical', ticks = yticks )  
        cbar.ax.set_yticklabels(xticks)
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        cbar.set_label('$Radius(kpc)$', fontsize=30 )
    
    ###########################################################################################################################
    ###########################################################################################################################    
    elif mode == 'each_in_5':
        
        fig = plt.figure(figsize=(15,10))
        #  theta and phi
        ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2, projection = 'polar')
        ax2 = plt.subplot2grid((2,6), (0,2), colspan=2, projection = 'polar')
        ax3 = plt.subplot2grid((2,6), (0,4), colspan=2, projection = 'polar')
        ax4 = plt.subplot2grid((2,6), (1,1), colspan=2, projection = 'polar')
        ax5 = plt.subplot2grid((2,6), (1,3), colspan=2, projection = 'polar')
        
        axs = [ax1,ax2,ax3,ax4,ax5]
        for i in range(5):
            

            axs[i].set_rmin(0)
            axs[i].set_rmax(120)
            axs[i].set_rlim(0,120)
            axs[i].set_rlabel_position(-22.5)
            #axs[i].set_origin(0)
            axs[i].plot([phi],[theta], marker = 'o', c = 'black') # plot disk
            axs[i].set_rticks([30,45,60,90,120])
            axs[i].grid(True)
            axs[i].plot([90],[0], marker = '*', c = 'y', markersize = 15) # plot minor axis
            axs[i].plot([0],[90], marker = 's', c = 'b', markersize = 15) # plot mid axis
            axs[i].plot([(np.pi/180.)*90],[90], marker = '^', c = 'g', markersize = 15) # plot major axis
            #axs[i].plot((np.pi/180.)*110,0, linewidth = 0)#
            axs[i].set_title(radii[i])
            #axs[i].margins(0)
            #axs[i].set_rlabel_position(95)
            #axs[i].legend()
    fig.suptitle(title, fontsize=16)        
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()    


# given axes from volumetric density, returns the estimated axes of the potential.
def get_potential_axes(axes):
    A = axes[0]/axes[2]
    B = axes[1]/axes[2]
    
    A = 3*(A-1)+1
    B = 3*(B-1)+1
    
    axes[1] = axes[0]*B/A
    axes[2] = axes[0]/A
    
    return axes        
    
    
# given axes from volumetric density, returns the estimated axes of the potential.
def get_potential_axes2(axes):
    A = axes[0]/axes[2]
    B = axes[1]/axes[2]
    
    A = A*np.sqrt(2*A**2-1)
    B = B*np.sqrt(2*B**2-1)
    
    axes[1] = axes[0]*B/A
    axes[2] = axes[0]/A
    
    return axes 
    
    
# given axes from volumetric density, returns the estimated axes of the potential.
def get_potential_axes3(axes):
    A = axes[1]/axes[0]
    B = axes[2]/axes[0]
    
    A = A*np.sqrt(2*A**2-1)
    B = B*np.sqrt(2*B**2-1)
    
    axes[1] = A*axes[0]
    axes[2] = B*axes[0]
    
    return axes        
        



main()
'''
quant = T_MHD - T_DM# Some measure of symmetry
    # np.abs(2*(0.5 - T_DM)) is 1 when object is axisymetric and 0 when it is maximally triaxial
    name = 'Delta Triaxiality'
    
    
    import matplotlib.pyplot as plt
    #plt.plot(axesMHD[:,0,1]/axesMHD[:,0,0],axesMHD[:,0,2]/axesMHD[:,0,0], 'ro' )
    #plt.plot(axesDM[:,0,1]/axesDM[:,0,0],axesDM[:,0,2]/axesDM[:,0,0], 'bo' )
    #plt.show()
    
    
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
    
    with PdfPages("./pics/alignments_level3"+".pdf") as pdf:
        halos = range(nhalos)
        #halos = [11,10]
        for j in halos:
            # Graphics
            rads = (axesMHD[j,:,0]*axesMHD[j,:,1]*axesMHD[j,:,2])**(1./3.)
            polar_graphs(pdf,vecsMHD[j],disk_vec[j],radii,"halo_"+str(j), radii_vals= rads, mode = 'all_in_1')
        
            
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
            
    
    
    # We study the axial ratios between different measurements
    #halos = range(30)
    halos = range(nhalos)
    ba_err= []
    ca_err= []
    for j in halos:
        # Text 
        print "______________________________________________________________"
        print "____________" + "halo_"+str(j)+"______________________________"
        print "______________________________________________________________\n"
        
        for i in range(5):
                print "############        " + radii[i] +"             ###############\n\n"
                print "Volume"
                print "b/a = " + str(axesMHD[j,i,1]/axesMHD[j,i,0]) 
                print "c/a = " + str(axesMHD[j,i,2]/axesMHD[j,i,0]) + "\n"
                
               
                print "IsoDensity"
                print "b/a = " + str(axesDen[j,i,1]/axesDen[j,i,0]) 
                print "c/a = " + str(axesDen[j,i,2]/axesDen[j,i,0]) + "\n"
                
                print "IsoPotential"
                print "b/a = " + str(axesPot[j,i,1]/axesPot[j,i,0]) 
                print "c/a = " + str(axesPot[j,i,2]/axesPot[j,i,0]) + "\n"
                
                
                print "IsoPotential Formula Isodensity"
                axesPot1 = np.copy(axesPot)
                axesPot1[j,i] = get_potential_axes(np.copy(axesPot[j,i]))
                print "b/a = " + str(axesPot1[j,i,1]/axesPot1[j,i,0]) 
                print "c/a = " + str(axesPot1[j,i,2]/axesPot1[j,i,0]) + "\n"
                
                ba_err.append(abs((axesPot1[j,i,1]/axesPot1[j,i,0])-(axesMHD[j,i,1]/axesMHD[j,i,0]))/(axesMHD[j,i,1]/axesMHD[j,i,0])) 
                ca_err.append(abs((axesPot1[j,i,2]/axesPot1[j,i,0])-(axesMHD[j,i,2]/axesMHD[j,i,0]))/(axesMHD[j,i,2]/axesMHD[j,i,0]))
                
                
                
                print "\nTriaxiality = " + str(T_MHD[j,i])

    ba_err = np.array(ba_err).reshape(nhalos,5)
    ca_err = np.array(ca_err).reshape(nhalos,5)                              
    print np.mean(ba_err,axis=0),np.std(ba_err,axis=0)
    print np.mean(ca_err,axis=0),np.std(ca_err,axis=0)    

'''


















