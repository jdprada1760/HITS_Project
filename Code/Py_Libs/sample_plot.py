import matplotlib
matplotlib.use('Agg')

from arepo import *
import numpy as np
import matplotlib.pyplot as plt
import yt
from yt.units import parsec, Msun

# Specification of the simulation 
path = "/hits/universe/GigaGalaxy/level5/halo16/output" 
snap = "/snapshot_063.hdf5"
subSnap = "/fof_subhalo_tab_063.hdf5"

# Loads the simulation to make some graphics
sn = Simulation(path+snap, parttype=[1])
# Loads groups to select the principal halo
subSn = Subfind(path+subSnap)

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

ppx,ppy,ppz = pos.T
# Creates data to work with in YT
data = {'particle_position_x': pos[:,0],
        'particle_position_y': pos[:,1],
        'particle_position_z': pos[:,2],
        'particle_mass': sn.MassTable[1]*np.ones(len(pos))}

bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])
ds = yt.load_particles(data, length_unit=1000*parsec, mass_unit=1e10*Msun, n_ref=64, bbox=bbox)
ad = ds.all_data()

# This is generated with "cloud-in-cell" interpolation.
cic_density = "deposit", "all_cic"
# Neares Neighbor
nn_density = "deposit", "all_density"
nn_deposited_mass = "deposit", "all_mass"
particle_count_per_cell = "deposit", "all_count"

# Choose a center for the render.
c = [0, 0, 0]



'''
# Create the off axis projection.
# Setting no_ghost to False speeds up the process, but makes a
# slighly lower quality image.

p = yt.ProjectionPlot(ds, 'z', nn_density, center=c, weight_field=nn_density, width=(2*rad, 'kpc'))
#p = yt.ParticlePlot(ds, 'particle_position_x','particle_position_y', weight_field=nn_density, width=(2*rad,2*rad))

# Title
p.annotate_title('Density Plot')
# Changes CMap
p.set_cmap(nn_density,'inferno')
# Draws density contour
p.annotate_contour(nn_density, ncont = 10, take_log = True )
# Plots center of the halo (potential) 
p.annotate_marker((0,0), coord_system='plot',plot_args={'color':'blue','s':500})

p.save("OtherAttempt.png")
#p.display()

##############################################################################################

# Gets matplotlib figure,axes
mpl = p.plots[nn_density]
mpl.axes.plot([0],[0], marker = 'o')
# Draws ellipse
from matplotlib.patches import Ellipse
elli = Ellipse(xy=[0,0],width = rad, height = 0.5*rad, fill = False, linewidth = 2.5)  
# Filter box
mpl.axes.add_artist(elli)
#mpl.figure.savefig("Otro.png")

#############################################################################################


mpl.figure.savefig("Alternative.png")

#image = yt.off_axis_projection(ds, c, L, W, Npixels, ('deposit','all_cic'), no_ghost=True)
#image.set_cmap(('deposit','all_cic'),'inferno')
'''






'''

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


Ni = None
# Generates ellipse
a,b,c = 4.,3.,2.
pos = gen_ellipse(a,b,c,100000)#-np.array([2,.1,1])
#np.array(np.random.random(3)*np.array([a,b,c])*0.5)
# Random rotation
theta = np.random.random()*np.pi 
Rot = np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
pos = Rot.dot(pos.T).T
theta = np.random.random()*np.pi 
Rot = np.array([[np.cos(theta),0,-np.sin(theta)],[0,1,0],[np.sin(theta),0,np.cos(theta)]])
pos = Rot.dot(pos.T).T
rad = 5
a,b,c,pos = recalculate_Semiaxes(rad,rad,rad,pos)
from matplotlib.patches import Ellipse
elli = Ellipse(xy=[0,0],width = 2*a, height = 2*b, fill = False)  
# Filter box
ind = np.where((abs(pos[:,0])< 1.5*a) & (abs(pos[:,1])< 1.5*a)& (abs(pos[:,2])< 1.5*a))[0]
npos = pos[ind,:]
# Plot
fig = plt.figure(0)
ax = fig.add_subplot(111, aspect='equal')
ax.add_artist(elli)
markersize = 0.1

if(Ni is None):
    alpha = 0.2
    markersize = 0.5
#alpha = (10.0**(-1.0*len(npos)/Ni))
alpha = 0.7
#abel = "r="+'{:.2e}'.format(float(rad))
ax.plot([0],[0], marker = 'o')
ax.plot(npos[:,0], npos[:,1], marker = '.', linewidth = 0, alpha = alpha, markersize = markersize)
ax.set_xlim(-1.5*a,1.5*a)
ax.set_ylim(-1.5*a,1.5*a)
title = "Teo:0.75_0.5_____Exp:"
title += '{:.4e}'.format(float(b/a))+"_"
title += '{:.4e}'.format(float(c/a))
#label += "\n c="+'{:.2e}'.format(float(c))
ax.set_title(title)
plt.savefig("../Plots/Test/Elli_5_Normal.png")
print(b/a,c/a,'|__|',3./4.,2./4.)
print(np.abs((b/a)-(3./4.))/(3./4.),np.abs((c/a)-(2./4.))/(2./4.))
'''
