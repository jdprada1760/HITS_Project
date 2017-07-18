import matplotlib
matplotlib.use('Agg')

from arepo import *
import numpy as np
import matplotlib.pyplot as plt

# Specification of the simulation 
path = "/hits/universe/GigaGalaxy/level5/halo16/output" 
snap = "/snapshot_063.hdf5"
subSnap = "/fof_subhalo_tab_063.hdf5"

# Loads groups to select the principal halo
subSn = Subfind(path+subSnap)
filterr = [Halo(subSn, halo = 0, subhalo = subSn.GroupFirstSub)]

# Loads the simulation to make some graphics
sn = Snapshot(path+snap, parttype=[1], filter=filterr, combineFiles=True, verbose=True)

# r the radius of the halo (R_mean200)
rad = subSn.group.Group_R_Mean200[0]/2

# Calculates the center of mass f the subhalo (initial)
CM = subSn.SubhaloCM[0]
# 
CM2 = subSn.GroupCM[0]-CM
# The particle positions
pos =  sn.pos
pos -= CM

# Initial number of particles
Ni = len(pos)

# Gets semiaxes given eigenvalues of Inertia and a radius
def calc_semiaxes(vals, rad):
    # Ellipse axial ratios
    q = np.sqrt(vals[1]/vals[0])
    s = np.sqrt(vals[2]/vals[0])
    # Calculates semiaxes of ellipsoid conserving volume
    a = ((rad**3.)/(q*s))**(1./3.)
    b = q*a
    c = s*a
    return a,b,c


# Recalculates semiaxes given the particle positions of the Halo and semiaxes before 
# Reference: https://arxiv.org/abs/1104.1566
# a,b,c are The semiaxes defining the elipsoid that encloses the considered particles
# cm is the center of mass to calculate Inertia
# returns the new semiaxes, cm and rotated particles along the semiaxes
def get_Semiaxes(a,b,c,pos):
    # Filter particles inside determined ellipsoid
    aux = (pos[:,0]/a)**2+(pos[:,1]/b)**2+(pos[:,2]/c)**2
    npos = pos[np.where(aux <= 1)[0],:]
    # Recalculates center of mass and centers particles (filtered and non-filtered)
    cm = npos.mean(axis = 0)
    npos = npos - cm
    pos = pos - cm
    # Axial ratios
    q = b/a
    s = c/a
    # The inertia tensor
    I = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            # The distance measure (usual distance when treating a circle)
            d = np.sum(npos[:,0]**2+(npos[:,1]/q)**2+(npos[:,2]/s)**2)
            I[i,j] = np.sum(npos[:,i]*npos[:,j])/d
    # Solves eigensystem and order values
    vals,vecs = np.linalg.eig(I)
    vecs = vecs[:,(-vals).argsort()]
    vals[::-1].sort()
    # Gets semiaxes
    a,b,c = calc_semiaxes(vals,rad)
    # Rotates positions to diagonal basis
    pos = (vecs.T).dot(pos.T).T
    return a,b,c,pos
    
# Calls get_semiaxes iteratively reshaping iteratively with previous result starting from initial calculation
# Initial semiaxes a,b,c
# n_it = number of iterations
# cm the initial center of mass
def recalculate_Semiaxes(a,b,c,pos,n_it=100):
    # Iterates over n_it
    for i in range(n_it):
        # Calculates inertia tensor
        a,b,c,pos = get_Semiaxes(a,b,c,pos)
    return a,b,c,pos

def plotHalo(filename,a,b,c,pos):
    rad = (a*b*c)**(1./3.)
    print(a,b,c,rad)
    # Draws ellipse
    from matplotlib.patches import Ellipse
    elli = Ellipse(xy=[0,0],width = a, height = b, fill = False)   
    # Filter box
    ind = np.where((abs(pos[:,0])< rad) & (abs(pos[:,1])< rad) & (abs(pos[:,2])< rad))[0]
    npos = pos[ind,:]
    # Plot
    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')
    ax.add_artist(elli)
    if(len(npos) == 0):
        print("No particles")
        alpha = 1
    else:
        alpha = (10.0**(-1.0*len(npos)/Ni))
    ax.plot(npos[:,0], npos[:,1], marker = '.', linewidth = 0, alpha = alpha, markersize = 0.1, label = "r="+'{:.2e}'.format(float(rad)))
    ax.legend()
    plt.savefig(filename)
    plt.close()
    
# Looking if it works
#recalculate_Inertia(pos,rad,50)

# Logspace for radius
xrad = np.logspace(np.log10(0.01*rad),np.log10(2*rad),num=10, endpoint = True)[::-1]
#print(0.01*rad,2*rad)
#print(xrad)
a,b,c = xrad[0],xrad[0],xrad[0]
for i in range(len(xrad[:-1])):
    rad = (a*b*c)**(1./3.)
    print(xrad[i],rad)
    a,b,c,pos = recalculate_Semiaxes(a,b,c,pos)
    plotHalo("lvl5_h16_rm200_"+str(i)+".png",a,b,c,pos)
    # Uses last semiaxes to simplify further calculations by rescaling
    scale = xrad[i+1]/xrad[i]
    a,b,c = scale*a,scale*b,scale*c
    
