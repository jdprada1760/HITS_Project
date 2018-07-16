#import matplotlib
#matplotlib.use('Agg')

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
    gas_den = par['MassGas']#/(GasDiskRad**2)     # Mass densities
    star_den = par['MassStars']#/(StarDiskRad**2)
    BH_den = par['MassBH']#/(StarDiskRad**2)
    bar_frac = (par['MassGas']+par['MassStars'])/(par['MassDM']+par['MassGas']+par['MassStars']) # The baryonic fraction
    
    '''
    import matplotlib.pyplot as plt
    
    plt.plot(StarDiskRad,bar_frac, 'ro')
    plt.xlabel('Star Rad')
    plt.ylabel('Baryonic Fraction')
    plt.savefig('./pics/Baryons_plots/StarRad_barfrac.png')
    plt.close()
    
    plt.plot(np.log10(star_den+gas_den),bar_frac, 'ro')
    plt.xlabel('Star mass')
    plt.ylabel('Baryonic Fraction')
    plt.savefig('./pics/Baryons_plots/StarMass_barfrac.png')
    plt.close()
    '''
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
    #y = (asph_MHD - 0*asph_DM)[:,1] # 0 -> only calculated at 0.125 Rvir
    
    #y1 = axesMHD[:,0,1]/axesMHD[:,0,0] # this is b/a MHD for all haloes
    y = ((axesMHD[:,0,2]/axesMHD[:,0,0]))**(-1) #- (axesDM[:,0,2]/axesDM[:,0,0])) # c/a MHD

    # Just to see if the increase is related to the unperturbed state of the halo
    x1 = axesDM[:,0,1]/axesDM[:,0,0] #b/a DM (q)
    x2 = axesDM[:,0,2]/axesDM[:,0,0] #c/a MHD (s)
    x1 = 1./x1
    x2 = 1./x2
    X = np.array([GasDiskRad,StarDiskRad,np.log10(gas_den),np.log10(star_den),np.log10(BH_den),bar_frac,x1,x2]).T
        
    Xnames = ['Gas Radius','StarDisk Radius','Gas density','Star density','BH density',
                'Baryonic fraction','q','s']
    Xlog = [0,0,1,1,1,0,0,0,0] # To know if scales are logarithm
    
    
    # Halos to discard
    #i_discard = np.where(y<0.1)[0]
    i_discard = [10]
    badX = X[i_discard]
    bady = y[i_discard]
    X = np.delete(X,i_discard,axis=0)
    y = np.delete(y,i_discard)
    
    #SVReg(X,(y),Xnames=Xnames, yname = 'Delta sphericity', Xlog = Xlog, badX=badX, bady=bady)
      
    # Mean values
    meanX = np.mean(X, axis = 0)
    # STD but with the signs of the coeficients
    sigmaX = np.std(X, axis = 0)
    
    #X = ((X-meanX)/sigmaX)


    MCMC(X,y,X_names=Xnames)

    
    
    '''
    print meanX[1]
     
    from scipy.optimize import curve_fit
    def fun( x,a,b,c,n): 
        return a/((x-b)**n)+c
    popt,pcov = curve_fit(fun, X[:,1],X[:,-1])
        
    #scatter3d(X[:,-2], X[:,-1], y, y)    
    plt.plot(X[:,1],X[:,-1],'ro')
    xspc = np.linspace(min(X[:,1])-0.1, max(X[:,1])+0.1, 100)
    print popt
    plt.plot(xspc,fun(xspc,popt[0],popt[1],popt[2],popt[3]))   

    plt.show()
    '''



def MCMC(X,y,X_names='',yname=''):
    
    # Basic imports
    import emcee
    
    
    
    ##########################
    ##### SETUP FUNCTIONS ####
    ##########################
    print X.shape
    X = X.T
    
    # The model to fit parameters
    # Asume linear model and power law for baryionic fraction
    # Encode arguments of this function in theta
    def model(theta,X):

        theta = np.array(theta)
        #print np.multiply(theta[:-1][:,np.newaxis],X)
        val = np.sum(np.multiply(theta[:-1][:,np.newaxis],X),axis=0)+ theta[-1]
        #print val
        return val

    # Define the ln of the likelihood.
    # Just assume exponential of the squared error
    def lnlike(theta, X, y, yerr=1):
        
        y_mod = model(theta,X)
        inv_sigma2 = 1.0/(yerr**2) #+ model**2*np.exp(2*lnf))      -> This is more sophisticated
        return -0.5*(np.sum((y-y_mod)**2*inv_sigma2 - np.log(inv_sigma2)))
    
    # The log of the prior probability
    # All variables are free to choose (mainly slopes)
    def lnprior(theta,X):
        
        #a0,a1,c = theta
        if False : # keep power term positive
            return -np.inf
        return 0 # no restrictions
    
    # The log of the final probability of the model    
    def lnprob(theta, X, y, yerr=1):
        lp = lnprior(theta,X)
        if not np.isfinite(lp): # Not needed but anyway jic
            return -np.inf
        return lp + lnlike(theta, X, y, yerr)
    
    # Now we find the maximum likelihood
    import scipy.optimize as op
    nll = lambda x,*t: model([t[i] for i in range(len(X)+1)],x)
    print len(X.T),len(y)
    x = np.transpose(X)
    popt,pcov = op.curve_fit(nll,X,y,p0=0*np.ones(len(X)+1))
    
    import matplotlib.pyplot as plt     
    #scatter3d(X[-2], X[-1], y, y)    
    yopt = model(popt,X)
    xlabs = X_names
    for i in range(len(X)):
        plt.plot(X[i],y,'bo')
        plt.plot(X[i],yopt,'ro')
        plt.ylabel(r'$ \frac{a}{c}$')
        plt.xlabel(xlabs[i])       
        plt.savefig('./pics/Baryons_plots/'+xlabs[i]+'.png')
        plt.close()
    #xspc = np.linspace(min(X[:,1])-0.1, max(X[:,1])+0.1, 100)
    #print popt
    #plt.plot(xspc,fun(xspc,popt[0],popt[1],popt[2],popt[3]))   

    
    
    
    
    
    print popt
    # We initialize some walkers around the maximum likelihood
    ndim, nwalkers = len(X)+1, 100
    pos = [popt + (1e-6)*(1-2*np.random.randn(ndim)) for i in range(nwalkers)]
    
    # Set up the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(X, y, 1.0))
    
    # Run nwalks 
    sampler.run_mcmc(pos, 1000)
    
    # Retrieve chains after stabilization
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    
    #print result
    import corner    
    fig = corner.corner(samples, labels= X_names+['$Intercept$'], truths = popt)
    
    results = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))
    
    for el in zip(results,X_names+['Intercept']):
        print el
    
    # Tries to create directory
    path = './pics/Baryons_plots/'
    try:
        os.makedirs(path)
    except OSError:
        print("Directory already exists")
    fig.savefig(path+'MCMC_mixed.png')
    
    

    

# Performs a Support vector regression according the given yparams and xparams
# Plots important projections over the x space to see accuracy of fit.
# Returns parameters of the fit.    
def SVReg(X,y,Xnames='',yname='',Xlog = None, badX=None, bady=None):        
    
    import matplotlib.pyplot as plt
    # Import scikit learn
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import Lasso

    # And random forest
    rf = RandomForestClassifier(n_estimators=100)
    # And LassoCV
    lasso = Lasso(alpha =0.01)
    
    
    # Mean values
    meanX = np.mean(X, axis = 0)
    # STD but with the signs of the coeficients
    sigmaX = np.std(X, axis = 0)
    
    X = ((X-meanX)/sigmaX)
    badX =  ((badX-meanX)/sigmaX)
    
    # Labels for random forest
    mid1 = np.percentile(y,33)
    mid2 = np.percentile(y,66)
    labels = (y>mid2).astype(np.int)
    #labels+= (y>mid1).astype(np.int)
    print(labels)
    
    # Asuming array x is of size (npoints,mparams) and array y is of size (npoints)

    lasso.fit(X,y)
    
    rf.fit(X,labels)
    
    fit_labs = rf.predict(X)
    
    # Linspaces for each variable
    linspcX = np.array([np.linspace(min(x)-0.01,max(x)+0.01,100) for x in X.T])
    #print Xs
    
    y_data_lasso = lasso.predict(X) 
    
    for i in range(len(Xnames)):
        # Names and slopes
        print Xnames[i],
        print '\n',
        print 'Slope',lasso.coef_[i]
        print 'Feature Importance',rf.feature_importances_[i]
        print 'Feature correlation',np.corrcoef(X[:,i],y)[0,1],'\n'
        
        #plt.plot(X[:,1],X[:,i], 'bo')
        
        # Plot dataset projection over variable 
        
        plt.plot(X[:,i],y,'bo')
        
        plt.plot(badX[:,i],bady,'ro')
        
        #plt.hlines(mid, min(linspcX[i]), max(linspcX[i]), linestyle='-.')
        # We produce a generic space of with the linspace of the projected variable
        # and the mean in the other variables.
        # We also take the uppper and lower bounds in the other variables
        
        #X_proj_mean = np.zeros((100,len(meanX)))
        #X_proj_mean[:,i] = linspcX[i]

        #yproj_svr = svr.predict(X_proj_mean)
        #yproj_rf = rf.predict(X_proj_mean)
        #yproj_lasso = lasso.predict(X_proj_mean)
        
        #plt.plot(linspcX[i],yproj_svr,'black',linewidth=2,label='SVR') # Mean
        #plt.plot(X[:,i],y_data_svr,'g*',label='SVR_data') # Mean
        
        
        
        #plt.plot(linspcX[i],yproj_lasso,'r',linewidth=2,label='Lasso') # Mean
        #plt.plot(X[:,i],y_data_lasso,'m^',linewidth=2,label='Lasso_data') # Mean

        plt.legend(loc=0)
        
        plt.xlabel(Xnames[i])
        #plt.show()
        plt.close()
    
    print "FIN"     


def scatter3d(x,y,z, cs, colorsMap='jet'):
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cmx
    from mpl_toolkits.mplot3d import Axes3D
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.show()
    
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
            
            #print "Minor: "
            #print vecs[0,[1,0,2],:].dot(vecs[i,2,:])
            #print to_Cartesians(mint,minp)
            #print "Mid: "
            #print vecs[0,[1,0,2],:].dot(vecs[i,1,:])
            #print to_Cartesians(midt,midp)
            #print "Major: "
            #print vecs[0,[1,0,2],:].dot(vecs[i,0,:])
            #print to_Cartesians(majt,majp)
            #print "____________________________________________________"
 
            
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


















