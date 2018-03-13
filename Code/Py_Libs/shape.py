import matplotlib
matplotlib.use('Agg')

from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import yt
from yt.units import parsec, Msun

DMass = 0


# Loads given snapshot of the simulation given the needed parameters
# lvl -> The level of the simulation
# halo -> The name of the halo
# snap -> Name of snap to load
# subSnap -> Name of the file where group and subgroup information is stored
# return pos -> The centered position of the main structure particles (filter substructures)
# return rad -> The virial radius 
def loadSim(lvl,halo,snapnum):
    global DMass

    # Path of the simulation
    path = "/hits/universe/GigaGalaxy/"+lvl+"/"+halo+"/output" 

    # Loads groups to select the principal halo
    subSn = Subfind(path,snapnum)

    # filterr = [Halo(subSn, halo = 0)]
    # filterr = [Halo(subSn, halo = 0, subhalo = subSn.GroupFirstSub[0])]

    # Loads the simulation
    #sn = Snapshot(path+snap, parttype=[1], filter=filterr, combineFiles=True, verbose=True)
    sn = Snapshot(path,snapnum, parttype=[1], combineFiles=True, verbose=True)
    
    # DM Mass
    DMass = sn.MassTable[1]
    
    # r the radius of the halo (R_mean200)*1000 to obtain it in kpc
    radiii = 1000.*subSn.group.Group_R_TopHat200[0]
    rad = 1000.*subSn.group.Group_R_Crit500[0]
    print("Rad_Crit_500 Vs Rad_Top_200: " + str(rad)+"   "+str(radiii))
    # Some galactic disk properties of the principal subhalo
    diskrads = 1000*subSn.SubhaloHalfmassRadType[0]
    
    # The redshift and scale factor
    redshift = sn.Redshift
    scalefactor = sn.Time

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

    # Creates a dictionary of important quantities about the studied
    dictio = {}
    dictio['Redshift'] = redshift
    dictio['Scale factor'] = scalefactor
    dictio['rad_gas'] = diskrads[0]
    dictio['rad_dm'] = diskrads[1]
    dictio['rad_stars'] = diskrads[4]
    dictio['rad_bh'] = diskrads[5]
  
    
    return pos,rad,dictio


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

# Plots Halo with YT
def plotHalo(filename,a,b,c,pos,lvl,halo):
    rad = (a*b*c)**(1./3.)
    # Filter box
    ind = np.where((abs(pos[:,0])< 1.5*a) & (abs(pos[:,1])< 1.5*a)& (abs(pos[:,2])< 2*a))[0]
    npos = pos[ind,:]
    ppx,ppy,ppz = npos.T

    # Creates data to work with in YT
    data = {'particle_position_x': ppx,
            'particle_position_y': ppy,
            'particle_position_z': ppz,
            'particle_mass': DMass*np.ones(len(npos))}

    bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])
    ds = yt.load_particles(data, length_unit=1000*parsec, mass_unit=1e10*Msun, n_ref=64, bbox=bbox)
    ad = ds.all_data()

    # This is generated with "cloud-in-cell" interpolation.
    #cic_density = "deposit", "all_cic"
    # Neares Neighbor
    nn_density = "deposit", "all_density"
    #nn_deposited_mass = "deposit", "all_mass"
    #particle_count_per_cell = "deposit", "all_count"

    # Choose a center for the render.
    c = [0, 0, 0]

    p = yt.ProjectionPlot(ds, 'z', nn_density, center=c, weight_field=nn_density, width=(3*a, 'kpc'))
    #p = yt.ParticlePlot(ds, 'particle_position_x','particle_position_y', weight_field=nn_density, width=(2*rad,2*rad))

    # Title
    p.annotate_title(filename.split('.')[0])
    # Changes CMap
    p.set_cmap(nn_density,'inferno')
    # Draws density contour
    p.annotate_contour(nn_density, ncont = 10, take_log = True )
    # Plots center of the halo (potential) 
    p.annotate_marker((0,0), coord_system='plot',plot_args={'color':'blue','s':500})

    p.save("tmp.png")
    #p.display()

    ##############################################################################################

    # Gets matplotlib figure,axes
    mpl = p.plots[nn_density]
    mpl.axes.plot([0],[0], marker = 'o')
    # Draws ellipse
    from matplotlib.patches import Ellipse
    elli = Ellipse(xy=[0,0],width = 2*a, height = 2*b, fill = False, linewidth = 2.5, label = "r="+'{:.2e}'.format(float(rad)))  
    # Filter box
    mpl.axes.add_artist(elli)
    leg  = mpl.axes.legend(facecolor='w')
    #text = leg.get_texts()
    #plt.setp(text, color = 'w')
    #mpl.figure.savefig("Otro.png")

    #############################################################################################


    mpl.figure.savefig("../Plots/"+lvl+"/"+halo+"/"+filename)

    plt.close()

def main():
    #lvl = 'level5'
    lvl = 'level4_MHD'
    #lvl = 'level5_Durham' 
    #halonums = [20]
    halonums = range(1,11)
    #halonums = ['halo16','halo16_MHD','halo24','halo24_MHD','halo28', 'halo6_MHD','halo9','halo9_MHD']
    for j in halonums:
        print("--------------------------------------------------------------------------")
        print(j)
        halo = 'halo_'+str(j)
        #halo = j
        #snap = "/snapdir_127/"
        #subSnap = "/groups_127/"
        snapnum = 127
        #snapnum = 63
        #snapnum = 255
        pos,rvir = loadSim(lvl,halo,snapnum)
        rvir
        print("Radius of the simulation:  "+str(rvir))
        # Logspace for radius
        logdelta = 3.0/50.0
        logr_ini = -1
        # Keeps semiaxes
        semmiaxes = []
        x_rad = 10.0**(logr_ini)
        a,b,c = x_rad,x_rad,x_rad 
        for i in range(1000):
            a,b,c,posx = recalculate_Semiaxes(a,b,c,pos)
            x_rad = 10.0**(logr_ini+(i+1)*logdelta)
            if a > 0:
                semmiaxes.append([a,b,c])
                # Uses last result to improve convergence
                b,c = (x_rad/a)*b,(x_rad/a)*c
                a = x_rad
                pos = posx
                m = 1
                if i%m == 0:
                    title = lvl+"_"+halo+"_"+str(i//m)+".png"
                    #plotHalo(title,a,b,c,posx,lvl,halo)
            else:
                a,b,c = x_rad,x_rad,x_rad

            rad = (abs(a*b*c))**(1./3.)
            #print("________________________________________________________________")
            print(x_rad,rad)
            if( x_rad > 2.0*rvir ):
                break
        semmiaxes = np.array(semmiaxes)
        np.savetxt("../Plots/"+lvl+"/"+halo+"/"+"abc_"+lvl+"_"+halo+".txt",semmiaxes, delimiter = ',')
    '''
        if(len(semmiaxes) != 0):
            a,b,c = semmiaxes.T
            yvals = np.array([b/a,c/a,c/b])
            ylabel = ['b/a','c/a','c/b']
            xvals = (a*b*c)**(1./3.)
            fig,axs = ratiosPlots(xvals, yvals, ylabel, rvir)  
            plt.savefig("../Plots/"+lvl+"/"+halo+"/"+"Axial_ratios_vs_R__"+lvl+"_"+halo+".png")
            plt.close()
    '''


