import matplotlib
matplotlib.use('Agg')

from arepo import *
import numpy as np
import matplotlib.pyplot as plt

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
    #print(cm)
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
            #d = np.sum(npos[:,0]**2+(npos[:,1]/q)**2+(npos[:,2]/s)**2)
            d = np.sum(npos[:,0]**2+(npos[:,1]/1)**2+(npos[:,2]/1)**2)
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
def recalculate_Semiaxes(a,b,c,pos,n_it=10):
    # Iterates over n_it
    for i in range(n_it):
        # Calculates inertia tensor
        a,b,c,pos = get_Semiaxes(a,b,c,pos)
        #print(a,b,c)
    return a,b,c,pos

def plotHalo(filename,a,b,c,pos,lvl,halo):
    rad = (a*b*c)**(1./3.)
    # Draws ellipse
    from matplotlib.patches import Ellipse
    elli = Ellipse(xy=[0,0],width = a, height = b, fill = False)  
    # Filter particles inside determined ellipsoid
    aux = (pos[:,0]/a)**2+(pos[:,1]/b)**2+(pos[:,2]/c)**2
    npos = pos[np.where(aux <= 1)[0],:]
    # Recalculates center of mass and centers particles (filtered and non-filtered)
    cm = npos.mean(axis = 0) 
    #print(cm)
    # Filter box
    ind = np.where((abs(pos[:,0])< rad) & (abs(pos[:,1])< rad)& (abs(pos[:,2])< rad))[0]
    npos = pos[ind,:]
    #npos = pos
    # Plot
    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')
    ax.add_artist(elli)
    markersize = 0.1
    '''
    if(Ni is None):
        alpha = 0.2
        markersize = 0.5
    '''
    alpha = (10.0**(-1.0*len(npos)/Ni))
    label = "r="+'{:.2e}'.format(float(rad))
    #label += "\n b/a="+'{:.2e}'.format(float(b/a))
    #label += "\n c/a="+'{:.2e}'.format(float(c/a))
    #label += "\n c="+'{:.2e}'.format(float(c))
    ax.plot(npos[:,0], npos[:,1], marker = '.', linewidth = 0, alpha = alpha, markersize = markersize, label = label)
    ax.set_title(filename.split('.')[0])
    ax.legend()
    plt.savefig("../Plots/"+lvl+"/"+halo+"/"+filename)
    plt.close()

# Generates random elliptical ditribution, gaussian in radius given semiaxes
# a,b,c semiaxes
def gen_ellipse(a,b,c, n_points):
    pos = []
    x,y,z = 0,0,0
    pos.append([x,y,z])
    while len(pos) < n_points :
        x_new = np.random.normal(loc = x, scale = a/4.)
        y_new = np.random.normal(loc = y, scale = b/4.)
        z_new = np.random.normal(loc = z, scale = c/4.)
        alpha = dist(x_new,y_new,z_new,a,b,c)/dist(x,y,z,a,b,c)
        if alpha >= 1:
            x,y,z = x_new,y_new,z_new
        else:
            if np.random.random()<alpha:
                x,y,z = x_new,y_new,z_new
        pos.append([x,y,z])
    return np.array(pos)


# Elliptical gaussian distribution
def dist(x,y,z,a,b,c):
    return (2*np.pi)**(-3./2.)*a*b*c*np.exp(-((x/a)**2+(y/b)**2+(z/c)**2))

'''
Ni = None
# Generates ellipse
a,b,c = 4.,3.,2.
pos = gen_ellipse(a,b,c,100000)-np.array([2,5,3])
#np.array(np.random.random(3)*np.array([a,b,c])*0.5)
# Random rotation
theta = np.random.random()*np.pi 
Rot = np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
pos = Rot.dot(pos.T).T
theta = np.random.random()*np.pi 
Rot = np.array([[np.cos(theta),0,-np.sin(theta)],[0,1,0],[np.sin(theta),0,np.cos(theta)]])
pos = Rot.dot(pos.T).T
rad = 5
a,b,c,pos = recalculate_Semiaxes(rad,rad,rad,pos,n_it = 50)
plotHalo("noCM_Ellipse_5_(2,5,3).png",a,b,c,pos)
print(b/a,c/a,'|__|',3./4.,2./4.)
print(np.abs((b/a)-(3./4.))/(3./4.),np.abs((c/a)-(2./4.))/(2./4.))

# Plot
plt.plot(pos[:,0],pos[:,1], marker = '.', linewidth = 0, alpha = 0.1, markersize =0.5)
plt.xlim(-4*a,4*a)
plt.ylim(-4*b,4*b)
#plt.zlim(-1.5*c,1.5*c)
plt.savefig("Processed_Ellipse.png")
plt.close()
'''


# Specification of the simulation
lvl = 'level5' 
halo = 'halo28'
path = "/hits/universe/GigaGalaxy/"+lvl+"/"+halo+"/output" 
snap = "/snapshot_063.hdf5"
subSnap = "/fof_subhalo_tab_063.hdf5"

# Loads groups to select the principal halo
subSn = Subfind(path+subSnap)
filterr = [Halo(subSn, halo = 0, subhalo = subSn.GroupFirstSub)]

# Loads the simulation to make some graphics
sn = Snapshot(path+snap, parttype=[1], filter=filterr, combineFiles=True, verbose=True)

# r the radius of the halo (R_mean200)
rad = subSn.group.Group_R_Mean200[0]

# Calculates the center of mass f the subhalo (initial)
CM = subSn.SubhaloCM[0]
# 
CM2 = subSn.GroupCM[0]-CM
# The particle positions
pos =  sn.pos
pos -= CM

# Initial number of particles
Ni = len(pos)

# Logspace for radius
xrad = np.logspace(np.log10(0.04*rad),np.log10(0.5*rad),num=50, endpoint = True)[::-1]
#print(0.01*rad,2*rad)
#print(xrad)
# Keeps semiaxes
semmiaxes = []
a,b,c = xrad[0],xrad[0],xrad[0]
for i in range(len(xrad[:-1])):
    rad = (a*b*c)**(1./3.)
    print("________________________________________________________________")
    print(xrad[i],rad)
    a,b,c,pos = recalculate_Semiaxes(a,b,c,pos,n_it = 20)
    semmiaxes.append([a,b,c])
    if i%2 == 0:
        title = lvl+"_"+halo+"_"+str(i)+".png"
        plotHalo(title,a,b,c,pos,lvl,halo)
    # Uses last semiaxes to simplify further calculations by rescaling
    scale = xrad[i+1]/xrad[i]
    a,b,c = scale*a,scale*b,scale*c

semmiaxes = np.array(semmiaxes)
fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
a,b,c = semmiaxes.T
print(b/a)
print("__________________________________________________________")
print(c/a)
print("__________________________________________________________")
print(c/b)
ax0.plot((xrad[:-1]), b/a)
ax1.plot((xrad[:-1]), c/a)
ax2.plot((xrad[:-1]), c/b)
ax0.set_xscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')
ax0.set_ylim(0,1)
ax1.set_ylim(0,1)
ax2.set_ylim(0,1)
ax0.set_ylabel("b/a")
ax1.set_ylabel("c/a")
ax2.set_ylabel("c/b")
ax2.set_xlabel("log(R)")
plt.savefig("../Plots/"+lvl+"/"+halo+"/"+"Axial_ratios_vs_R__"+lvl+"_"+halo+".png")

