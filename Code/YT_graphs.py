import matplotlib
matplotlib.use('Agg')

from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import yt
from yt.units import parsec, Msun


# Plots Halo with YT
def plotHalo(filename,pos,lvl,halo,redshift):
    
    # Font size
    MEDIUM_SIZE = 20
    SMALL_SIZE = 20
    SSSMALL_SIZE = 11
    plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

    
    #rad = (a*b*c)**(1./3.)
    # Filter box
    #ind = np.where((abs(pos[:,0])< 1.5*a) & (abs(pos[:,1])< 1.5*a)& (abs(pos[:,2])< 2*a))[0]
    npos = pos[:,:]
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

    p = yt.ProjectionPlot(ds, 'z', nn_density, center=c, weight_field=nn_density, width=(0.2*(npos.max()-npos.min()), 'kpc'))
    #p = yt.ParticlePlot(ds, 'particle_position_x','particle_position_y', weight_field=nn_density, width=(2*rad,2*rad))

    # Title
    p.annotate_title(filename.split('.')[0])
    # Changes CMap
    p.set_cmap(nn_density,'inferno')
    # Draws density contour
    p.annotate_contour(nn_density, ncont = 10, take_log = True )
    # Plots center of the halo (potential) 
    p.annotate_marker((0,0), coord_system='plot',plot_args={'color':'blue','s':500, 'label' : "z="+'{:.2e}'.format(float(redshift))})
    p.set_figure_size(6)

    #p.save("tmp.png")
    #p.display()


    ##############################################################################################

    # Gets matplotlib figure,axes
    mpl = p.plots[nn_density]
    
    #mpl.axes.plot([0],[0], marker = 'o', label = "z="+'{:.2e}'.format(float(redshift)))
    # Draws ellipse
    from matplotlib.patches import Ellipse
    elli = Ellipse(xy=[0,0],width = 2*0, height = 2*0, fill = False, linewidth = 0)  
    # Filter box
    mpl.axes.add_artist(elli)
    leg  = mpl.axes.legend(facecolor='w',loc = 0)
    #text = leg.get_texts()
    #plt.setp(text, color = 'w')
    #mpl.figure.savefig("Otro.png")

    #############################################################################################


    mpl.figure.savefig("../Plots/"+lvl+"/"+halo+"/"+filename)

    plt.close()

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
    pos =  1000.*sn.pos[:cut]
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




#lvl = 'level5'
lvl = 'level4_MHD'
#lvl = 'level5_Durham' 
#halonums = [20]
halonums = [16,21]
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
    
 
    for i in range(5):
                pos,rvir,info = loadSim(lvl,halo,snapnum)
                redshift = info['Redshift']
                title = str(i)+".png"
                plotHalo(title,pos,lvl,halo,redshift)
                snapnum -= 10
        
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

















