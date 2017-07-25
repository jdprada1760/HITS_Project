import matplotlib
matplotlib.use('Agg')

from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# Loads given snapshot of the simulation given the needed parameters
# lvl -> The level of the simulation
# halo -> The name of the halo
# snap -> Name of snap to load
# subSnap -> Name of the file where group and subgroup information is stored
# return pos -> The centered position of the main structure particles (filter substructures)
# return rad -> The virial radius 
def loadSim(lvl,halo,snap,subSnap):
    # Path of the simulation
    path = "/hits/universe/GigaGalaxy/"+lvl+"/"+halo+"/output" 

    # Loads groups to select the principal halo
    subSn = Subfind(path+subSnap)

    # filterr = [Halo(subSn, halo = 0)]
    # filterr = [Halo(subSn, halo = 0, subhalo = subSn.GroupFirstSub[0])]

    # Loads the simulation
    #sn = Snapshot(path+snap, parttype=[1], filter=filterr, combineFiles=True, verbose=True)
    sn = Snapshot(path+snap, parttype=[1], combineFiles=True, verbose=True)

    # r the radius of the halo (R_mean200)*1000 to obtain it in kpc
    rad = 1000.*subSn.group.Group_R_TopHat200[0]

    # Calculates the center of the halo (initial)*1000 to get it in kpc__
    # Center can be the center of mass or the minimum of potential
    # CM = 1000.*subSn.SubhaloCM[0]
    CM = 1000*subSn.group.GroupPos[0]

    # Gets the length of the main structure particles /// 0 - Principal halo, 1 - DM
    cut = subSn.group.GroupLenType[0,1]
    subCut = subSn.subhalo.SubhaloLenType[0,1]
    print("Group Len is: "+str(cut))
    print("SubGroup Len is: "+str(subCut))

    # The particle positions *1000 -> kpc
    pos =  1000.*sn.pos[:subCut]
    pos -= CM

    return pos,rad


# Gets semiaxes given eigenvalues of Inertia and a radius
def calc_semiaxes(vals, rad):
    # Ellipse axial ratios
    q = np.sqrt(vals[1]/vals[0])
    s = np.sqrt(vals[2]/vals[0])
    # Calculates semiaxes of ellipsoid conserving volume
    #a = ((rad**3.)/(q*s))**(1./3.)
    # Calculates semiaxes of ellipsoid keeping biggest eigenvalue equal to the radius
    a = rad    
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
    if len(npos) < 3000:
        return -1.,-1.,-1.,0.
    # Recalculates center of mass and centers particles (filtered and non-filtered)
    # cm = npos.mean(axis = 0)
    #print(cm)
    #npos = npos - cm
    #pos = pos - cm
    
    # Axial ratios
    q = b/a
    s = c/a
    # The inertia tensor
    I = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            # The distance measure (usual distance when treating a circle)
            d = (npos[:,0]**2+(npos[:,1]/q)**2+(npos[:,2]/s)**2)
            if len(np.where(d==0)[0])!=0:
                npos[np.where(d==0)[0]] = 1e-17*np.array([1,1,1])
                d = (npos[:,0]**2+(npos[:,1]/q)**2+(npos[:,2]/s)**2)
            #d = np.sum(npos[:,0]**2+(npos[:,1]/1)**2+(npos[:,2]/1)**2)
            I[i,j] = np.sum(npos[:,i]*npos[:,j]/d)
    # Solves eigensystem and order values
    vals,vecs = np.linalg.eig(I)
    ind = (-vals).argsort()
    vecs = vecs[:,ind]
    vals = vals[ind]
    # Gets semiaxes
    a,b,c = calc_semiaxes(vals,a)
    # Rotates positions to diagonal basis
    pos = (vecs.T).dot(pos.T).T
    #print(vecs.T.dot(I).dot(vecs))
    return a,b,c,pos
    
# Calls get_semiaxes iteratively reshaping iteratively with previous result starting from initial calculation
# Initial semiaxes a,b,c
# epsilon -> Convergence tolerance
def recalculate_Semiaxes(a,b,c,pos,epsilon=1e-6):
    # Iterates over n_it
    for i in range(1000):
        # Calculates inertia tensor
        an,bn,cn,pos = get_Semiaxes(a,b,c,pos)
        # Stops if something bad happens
        if an < 0:
            return -1.,-1.,-1.,0.
        # ratio between semiaxes in consecutive iterations
        ratio = abs((((an-a)/a +(bn-b)/b +(cn-c)/c))/3.0)
        a,b,c = an,bn,cn
        # Stops if convergence is achieved within 1e-6
        if ratio < epsilon:
            print("Convergence achieved at: ",i)
            return a,b,c,pos
        #print(ratio,a,b,c)
    return a,b,c,pos

def plotHalo(filename,a,b,c,pos,lvl,halo):
    rad = (a*b*c)**(1./3.)
    # Draws ellipse
    from matplotlib.patches import Ellipse
    elli = Ellipse(xy=[0,0],width = 2*a, height = 2*b, fill = False)  
    # Filter box
    ind = np.where((abs(pos[:,0])< 1.5*a) & (abs(pos[:,1])< 1.5*a)& (abs(pos[:,2])< 1.5*a))[0]
    npos = pos[ind,:]
    # Plot
    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')
    ax.add_artist(elli)
    markersize = 0.15*(rvir/rad)
    '''
    if(Ni is None):
        alpha = 0.2
        markersize = 0.5
    '''
    #alpha = (10.0**(-1.0*len(npos)/Ni))
    alpha = 0.1
    label = "r="+'{:.2e}'.format(float(rad))
    #label += "\n b/a="+'{:.2e}'.format(float(b/a))
    #label += "\n c/a="+'{:.2e}'.format(float(c/a))
    #label += "\n c="+'{:.2e}'.format(float(c))
    ax.plot([0],[0], marker = 'o')
    ax.plot(npos[:,0], npos[:,1], marker = '.', linewidth = 0, alpha = alpha, markersize = markersize, label = label)
    ax.set_xlim(-1.5*a,1.5*a)
    ax.set_ylim(-1.5*a,1.5*a)
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


# Plots many graphics of ratios (must be between 0 and 1)
# yvals -> Fields to plot on x axis (may be more than one array)
# xvals -> Only one array for x axis
def ratiosPlots(xvals, yvals, ylabel, rvir):  
     fig, axs = plt.subplots(nrows=len(yvals))
     for ax,yval,ylab in zip(axs,yvals,ylabel):
         ax.plot(xvals,yval, marker = '.')
         ax.plot([rvir,rvir],[0,1])
         ax.set_xscale('log')
         # Plotting ratios
         ax.set_ylim(0,1)
         # Valid for all Milkyway-like galaxies
         ax.set_xlim(0.1,350)
         ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
         ax.set_ylabel(ylab)
     axs[-1].set_xlabel("log(R(kpc/h))")
     return fig,axs

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
'''
lvl = 'level5_Durham'
#lvl = 'level5' 
for j in range(25,26):
    print("-------------------------------------------------")
    print(j)
    halo = 'halo_'+str(j)
    snap = "/snapshot_255.hdf5"
    subSnap = "/fof_subhalo_tab_255.hdf5"
    pos,rvir = loadSim(lvl,halo,snap,subSnap)
    rvir
    print("Radius of the simulation:  "+str(rvir))
    # Logspace for radius
    logdelta = 3.0/50.0
    logr_ini = -1
    # Keeps semiaxes
    semmiaxes = []
    for i in range(1000):
        x_rad = 10.0**(logr_ini+i*logdelta)
        #x_rad =340.0
        a,b,c,posx = recalculate_Semiaxes(x_rad,x_rad,x_rad,pos)
        if a > 0:
            semmiaxes.append([a,b,c])
            if i%10 == 0:
                title = lvl+"_"+halo+"_"+str(i//10)+".png"
                plotHalo(title,a,b,c,posx,lvl,halo)
        rad = (abs(a*b*c))**(1./3.)
        #print("________________________________________________________________")
        print(x_rad,rad)
        if( x_rad > 2.0*rvir ):
            break
        # Uses last semiaxes to simplify further calculations by rescaling
        #scale = xrad[i+1]/xrad[i]
        #a,b,c = scale*a,scale*b,scale*c
    semmiaxes = np.array(semmiaxes)
    np.savetxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt",semmiaxes, delimiter = ',')
    if(len(semmiaxes) != 0):
        a,b,c = semmiaxes.T
        yvals = np.array([b/a,c/a,c/b])
        ylabel = ['b/a','c/a','c/b']
        xvals = (a*b*c)**(1./3.)
        fig,axs = ratiosPlots(xvals, yvals, ylabel, rvir)  
        plt.savefig("../Plots/"+lvl+"/"+halo+"/"+"Axial_ratios_vs_R__"+lvl+"_"+halo+".png")
        plt.close()
'''
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
lvl = 'level5_Durham'
with PdfPages('../Plots/'+lvl+'/AxialRatios.pdf') as pdf:
    for j in range(1,31):
        halo = 'halo_'+str(j)
        print(halo)
        path = "../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt"
        a,b,c = np.loadtxt(path,delimiter = ',').T
        snap = "/snapshot_255.hdf5"
        subSnap = "/fof_subhalo_tab_255.hdf5"
        pos,rvir = loadSim(lvl,halo,snap,subSnap)
        yvals = np.array([b/a,c/a,(a**2-b**2)/(a**2-c**2)])
        ylabel = ['b/a','c/a','T']
        xvals = (a*b*c)**(1./3.)
        fig,axs = ratiosPlots(xvals, yvals, ylabel, rvir)
        axs[0].set_title(halo)
        pdf.savefig(fig)
        plt.close()
