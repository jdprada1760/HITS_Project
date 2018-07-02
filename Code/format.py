########################################################################################################
##### THIS SCRITP IS FOR CALCULATING A SUMARY OF IMPORTANT QUANTITIES FOR DM HALO SHAPE             ####
##### GENERATES:                                                                                    ####
#####               PARAMETERS .CSV (GAS,STELLAR,BH MASSES, SCALE RADII & OTHER IMPORTANT PARAMS)   ####
#####               AXES & VECS.CSV (TRIAXIAL CHARACTERIZATION AT 4to5 DIFFERENT RADII)             ####
#####                               MEAN (VOLUMETRIC) DENSITY CHARACTERIZATION                      ####
#####                               ISOPOTENTIAL SHELL CHARACTERIZATION                             ####
#####                               ISODENSITY SHELL CHARACTERIZATION                               ####
#####                                                                                               ####
##### ALL THESE QUANTITIES ARE CALCULATED FOR DM AND MHD RUNS (WHEN POSSIBLE)                       ####
########################################################################################################



import numpy as np
from arepo import *
import gc
import csv

# Imports of C libraries
import ctypes as ct
inertia = ct.cdll.LoadLibrary('./C_Libs/inertia.so')

def main():
    '''
    lvl = 'level3_'
    snapnum = 63 # MHD snapshot
    snapnumDM = 127 # DM snapshot
    halo_nums = [6,16,21,23,24,27]

    '''
    lvl = 'level4_'
    snapnum = 127 # MHD snapshot
    snapnumDM = 127 # DM snapshot
    halo_nums = range(1,31)
    

    # The corresponding list and their filenames
    params = []
    fn_params = lvl+'params'+'.csv'
    axes_sum = []
    fn_axes = lvl+'axes'+'.csv'
    vecs_sum = []
    fn_vecs = lvl+'vecs'+'.csv'
    pot_axes = []
    fn_potaxes = lvl+'potential_axes'+'.csv'
    pot_vecs = []
    fn_potvecs = lvl+'potential_vecs'+'.csv'
    den_axes = []
    fn_denaxes = lvl+'density_axes'+'.csv'
    den_vecs = []
    fn_denvecs = lvl+'density_vecs'+'.csv'

    it = 0 # keeps track of iteration to avoid copying column names
    for i in halo_nums:
        
        halo = 'halo_'+str(i)
        
        print("--------------------------------------------------------------------")
        print(halo)
        print("--------------------------------------------------------------------")
        
        #####################################################################################
        # Path of the simulation
        path = "/hits/universe/GigaGalaxy/"+lvl+'MHD'+"/"+halo+"/output" 
        pathDM = "/hits/universe/GigaGalaxy/"+lvl+'DM'+"/"+halo+"/output" 

        # Loads groups to select the principal halo
        subSn = Subfind(path,snapnum) # Same MHD (star and DM particles)
        subSnDM = Subfind(pathDM,snapnumDM) # For the DMonly simulations

        # Loads simulation particcles (Stars)
        sn = Snapshot(path,snapnum, parttype=[4], combineFiles=True, verbose=False)# Stars 
        snH = Snapshot(path,snapnum, parttype=[1], combineFiles=True, verbose=False)# DM in MHD
        snDM = Snapshot(pathDM,snapnumDM, parttype=[1], combineFiles=True, verbose=False)# DM in DMonly
        
        #####################################################################################
        
        # Calculates the center of the halo (initial)*1000 to get it in kpc__
        # Center can be the center of mass or the minimum of potential
        CM = 1000*subSn.group.GroupPos[0]
        CM_DM = 1000*subSnDM.group.GroupPos[0]

        # Gets the length of the main structure particles /// 0 - Principal halo, 4 - Stars
        subCut = subSn.subhalo.SubhaloLenType[0,4]
        # Gets the length of the main structure particles /// 0 - Principal halo, 1 - DM
        subCutH = subSn.subhalo.SubhaloLenType[0,1]#MHD
        subCutDM = subSnDM.subhalo.SubhaloLenType[0,1]#DMonly
        

        # The particle positions *1000 -> kpc
        pos =  1000.*sn.pos[:subCut] # StarParticles
        posH = 1000.*snH.pos[:subCutH] # DM particles in MHD
        posDM = 1000.*snDM.pos[:subCutDM] # DM particles in DMonly
        
        # The potential!
        pot = snH.pot[:subCutH]
        
        print(len(pos),len(posH),len(posDM))
        
        pos -= CM
        posH -= CM
        posDM -= CM_DM
        
        pos = np.array(pos,dtype = np.float) # To make things work with the C library
        posH = np.array(posH,dtype = np.float)
        posDM = np.array(posDM,dtype = np.float)

        # The particle velocities and masses
        vel = sn.vel[:subCut]
        mass = sn.mass[:subCut]

        # The formation time (negative for wind particles)
        time = sn.GFM_StellarFormationTime[:subCut]
        
        # Radii: Rad crit 500 and Disk radius
        Rad500 = 1000.*subSn.group.Group_R_Crit500[0]
        RadDisk = 1000.*subSn.SubhaloHalfmassRadType[0,0]
        
        #print(Rad500,RadDisk)
        
        # We get only stars
        pos = pos[time>0]
        vel = vel[time>0]
        mass = mass[time>0]
        time = time[time>0]

        # We get only the 10% of youngest stars
        cut10 = np.percentile(time,10)
        pos = pos[time<cut10]
        vel = vel[time<cut10]
        mass = mass[time<cut10]
        time = time[time<cut10]

        # Calculates disk axis
        axis = get_axis(mass,pos,vel)
        
        ####################################################################################################
        ####################################################################################################
        # Creates a dictionary of important quantities about the studied
        vals = []
        keys = []
        
        vals.append(halo)
        keys.append('HaloID')
        vals.append(subSn.SubhaloSFR[0])
        keys.append('SFR')
        vals.append(1000*subSn.SubhaloHalfmassRadType[0,0]) # Gas half mass rad as disk size
        keys.append('GasDiskRad')
        vals.append(1000*subSn.SubhaloHalfmassRadType[0,4]) # Stellar
        keys.append('StellarDiskRad')
        
        # All quantities are calculated within the disk size
        vals.append(subSn.SubhaloMassInHalfRadType[0,0])
        keys.append('MassGas')
        vals.append(subSn.SubhaloMassInHalfRadType[0,1])
        keys.append('MassDM')
        vals.append(subSn.SubhaloMassInHalfRadType[0,4])
        keys.append('MassStars')
        vals.append(subSn.SubhaloMassInHalfRadType[0,5])
        keys.append('MassBH')
        vals.append(1000.*subSn.group.Group_R_Crit500[0])
        keys.append('HaloRadCrit500')
        
        vals = vals + list(axis)
        keys.append('vDisk1')
        keys.append('vDisk2')
        keys.append('vDisk3')
        
        # Cleans variables of simulation
        sn = 0
        snH = 0
        snDM = 0
        subSn = 0
        subSnDM = 0
        
        ####################################################################################################
        ####################################################################################################
        vals_axes = []
        vals_axesDM = []
        vals_axesPot = []
        vals_axesDen = []
        keys_axes = []
        
        vals_vecs = []
        vals_vecsDM = []
        vals_vecsPot = []
        vals_vecsDen = []
        keys_vecs = []    
        
        # We want the triaxial parameters at four radii
        
        # Radii at which we are going to measure axial ratios
        rads = [0.125*Rad500,0.25*Rad500,0.5*Rad500,Rad500,RadDisk]
        #print rads
        #print("mean", np.mean(posH))
        # Names to differentiate columns
        names = ['12','25','50','1','Disk']
         
        for radi,nami in zip(rads,names):
            print("Getting Volume Triaxiality MHD at "+str(radi))
            vals_axes,vals_vecs,keys_axes,keys_vecs=append_triaxial_params(vals_axes,vals_vecs,keys_axes,keys_vecs,posH,radi,nami) # For normal galaxies
            #print(vals_axes[-3:])
            print("Getting Volume Triaxiality DM at "+str(radi))
            vals_axesDM,vals_vecsDM,mockvar1,mockvar2=append_triaxial_params(vals_axesDM,vals_vecsDM,[],[],posDM,radi,nami)# For halos DM (just for comparison)
            #print(vals_axesDM[-3:])
            print("Getting IsoPotential Triaxiality MHD at "+str(radi))
            vals_axesPot,vals_vecsPot,mockvar1,mockvar2=append_triaxial_params(vals_axesPot,vals_vecsPot,[],[],posH,radi,nami,'isopotential',pot) # For MHD potential (for comparison)
            #print(vals_axesPot[-3:])
            print("Getting IsoDensity Triaxiality MHD at "+str(radi))
            vals_axesDen,vals_vecsDen,mockvar1,mockvar2=append_triaxial_params(vals_axesDen,vals_vecsDen,[],[],posH,radi,nami,'isodensity') # For MHD potential (for comparison)
            #print(vals_axesPot[-3:])
        ####################################################################################################
        ####################################################################################################
        
        # Fills files lists
        if it == 0:
            params.append(keys)
            axes_sum.append(keys_axes)
            pot_axes.append(keys_axes)
            den_axes.append(keys_axes)
            
            vecs_sum.append(keys_vecs)
            pot_vecs.append(keys_vecs)
            den_vecs.append(keys_vecs)
            
            
        params.append(vals)
        axes_sum.append(vals_axes)
        axes_sum.append(vals_axesDM)
        vecs_sum.append(vals_vecs)
        vecs_sum.append(vals_vecsDM)
        pot_axes.append(vals_axesPot)
        pot_vecs.append(vals_vecsPot)
        den_axes.append(vals_axesDen)
        den_vecs.append(vals_vecsDen)
        
        it += 1
        
        
        gc.collect() # Garbage collector does not clean vars from simulation loading
        

    # Save the corresponding lists
    save_csv(params,fn_params)
    save_csv(axes_sum,fn_axes)
    save_csv(vecs_sum,fn_vecs)
    save_csv(pot_axes,fn_potaxes)
    save_csv(pot_vecs,fn_potvecs)
    save_csv(den_axes,fn_denaxes)
    save_csv(den_vecs,fn_denvecs)





####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################    
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

# Saves list given by parameters in a "path/filename" also stated

def save_csv(list_save,filename):

    with open( '../Plots/'+filename, 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_NONE, quotechar='')
        wr.writerows(list_save)


# Appends to the list of values and keys the corresponding parameters of triaxiality
# given the differentiating string -name-
def append_triaxial_params(vals_axes,vals_vecs,keys_axes,keys_vecs,posH,rad,name, method='volume', potential = None):
    
    # Initialize variables
    axes = np.ones(3,dtype = np.float)*rad
    vecs = np.identity(3,dtype = np.float)
    
    # Calculation
    if method == 'volume':
        axes,vecs = get_ratios_vecs(posH,rad)
    elif method == 'isopotential':
        axes,vecs = triax_methods(posH,rad,method,potential)
    elif method == 'isodensity':     
        axes,vecs = triax_methods(posH,rad,method)
    
    
    
    vals_axes = vals_axes + list(np.copy(axes))
    # name formatting
    ax_names = ['a','b','c']
    keys_axes = keys_axes + [col+name for col in ax_names]
    
    vals_vecs = vals_vecs + list(np.copy(vecs.flatten()))
    # Name formatting
    vec_names = ['vax','vay','vaz','vbx','vby','vbz','vcx','vcy','vcz']    
    keys_vecs = keys_vecs + [col+name for col in vec_names]
    
    return vals_axes,vals_vecs,keys_axes,keys_vecs


# Calculates the angular momentum of the stars defined by the arrays as parameters
# Returns the normalized vector of that angular momentum
def get_axis(mass,pos,vel):
    L = np.array([0.0,0.0,0.0])
    L[0] = np.sum(mass*(pos[:,1]*vel[:,2]-pos[:,2]*vel[:,1]))
    L[1] = np.sum(mass*(pos[:,2]*vel[:,0]-pos[:,0]*vel[:,2]))
    L[2] = np.sum(mass*(pos[:,0]*vel[:,1]-pos[:,1]*vel[:,0]))
    
    return L/np.linalg.norm(L)


# Calculates the axial ratios and vectors of the DM halo at certain radius
# returns [a,b,c],[vec1,vec2,vec3]
def get_ratios_vecs(posH,rad):
    # Initialize variables
    axes = np.ones(3,dtype = np.float)*rad
    vecs = np.identity(3,dtype = np.float)
    
    # Ensures that we obtain an output radius near the corresponding radius
    rad_temp = 0
    l = 0
    while True:

        inertia.get_shape(ct.c_void_p(posH.ctypes.data),ct.c_int(len(posH)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))

        rad_temp = (abs(axes[0]*axes[1]*axes[2]))**(1./3.)
        
        # Advance
        if axes[0] <= 0:
            break
        elif abs(rad_temp-rad)/rad > 0.05 :
            axes = (axes[0]+.5*(rad-rad_temp))*axes/axes[0]
            #axes = np.ones(3,dtype = np.float)*(axes[0]+1.*(x_rad-rad))
            #vecs = np.identity(3,dtype = np.float)  
        else:
            break
            
        l += 1   
        #print("+++"+str(l)+"+++",rad,x_rad,axes[0])  
        if l > 30 :
            break
        
    return axes,vecs

# Given particles from posH and a radius rad and the respective potential of the particles
# Calculates the average potential (and std) at a certain radius \pm 10% 
# Gets particles that are within the mean \pm std and calculates triaxiality until convergence is achieved
# returns axes,vecs

def triax_methods(pos,rad, method = 'isodensity',pot = None):
    
    # Inertia tensor of this subset of particles
    # Initialize variables
    axes = np.ones(3,dtype = np.float)*rad
    vecs = np.identity(3,dtype = np.float)
    ratio = 10
    
    if method == 'isodensity':
        epsilon = 0.1
        tolerance = 1e-4
    elif method == 'isopotential':
        epsilon = 0.025
        tolerance = 1e-4
    
    for i in range(50):
        
        # First filter particles within a spherical shell to calculate potential
        d = np.linalg.norm((vecs.T).dot(pos.T).T*(axes/axes[0]),axis=1)/rad
        indices = np.where((d<1.0+epsilon)&(d>1-epsilon))[0] # Shell with 10%width of rad

        if method == 'isopotential':
            # Now we calculate the mean and desvest to capture only particles within these values
            mean_pot = np.mean(pot[indices])
            std_pot = np.std(pot[indices])
            
            #print("Mean:  "+str(mean_pot))
            #print("STD:  "+str(std_pot))
            # We filter particles with potential within mean \pm std
            indices = np.where((pot>(mean_pot-std_pot))&(pot<(mean_pot+std_pot)))[0] # We recycle this variable        
        

        # We calculate the eigensystem of the inertia tensor
        naxes,vecs = get_Eigensys(pos[indices], axes, vecs) # These axes are not normalized
        ratio = np.mean(abs((axes-rad*naxes)/(rad*naxes)))
        if ratio < tolerance:
            break
 
        # get normalized axes    
        axes = rad*naxes
    print("It:  " +str(i),"_____________________________ Conv: "+str(ratio))
    #print((axes[0]*axes[1]*axes[2])**(0.333), rad)    
 
        
    return axes,vecs
    


# calculates inertia tensor of particles given in pos taking into account the q(y->y/q) and s(z->z/s))
# Returns eigensys
def get_Eigensys(pos, axes, vecs):

    # Axial ratios
    q = (axes[1]/axes[0])
    s = (axes[2]/axes[0])
    
    # Rotates positions to diagonal basis
    
    # The inertia tensor
    I = np.zeros((3,3))
    
    # The distance measure (usual distance when treating a circle)
    d = np.linalg.norm((vecs.T).dot(pos.T).T*(axes/axes[0]),axis=1)**2
    zeros = np.where(d==0)[0]
    if len(zeros)!=0:
        #sp_pos[zeros] = 1e-17*np.array([1,1,1])
        d[zeros] = 1e17
        
    for i in range(3):
        for j in range(i,3):
            I[i,j] = np.sum(pos[:,i]*pos[:,j]/d)
            I[j,i] = I[i,j]
            
    # Solves eigensystem and order values
    vals,vecs = np.linalg.eig(I)
    ind = (-vals).argsort()
    vecs = vecs[:,ind]
    vals = vals[ind]
    # Gets semiaxes
    vals = np.sqrt(vals/vals[0])

    
    return vals,vecs
    
    
    
    
main()













