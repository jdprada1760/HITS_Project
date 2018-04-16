import matplotlib
matplotlib.use('Agg')

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

'''
#import yt
#from yt.units import parsec, Msun
from Py_Libs.shape import *

# C imports
import ctypes as ct
inertia = ct.cdll.LoadLibrary('./C_Libs/inertia.so')

# Simiulation specs
lvl = 'level3_MHD'
#lvl = 'level3_DM'
#halonums = range(1,31)
halonums = [6,16,21,23,24,27]
#snapnum = 127
snapnum = 63
#snapnum = 255

# Keeps semiaxes and eigenvecs
semmiaxes = []
eigenvecs = []
for j in halonums:
    # Loads simulation
    print("--------------------------------------------------------------------------")
    print(j)
    halo = 'halo_'+str(j)
    pos,rvir = loadSim(lvl,halo,snapnum)
    pos = np.array(pos,dtype = np.float)
    print("Radius of the simulation:  "+str(rvir))
    # Pointers as parameters for C_function
    axes = np.ones(3,dtype = np.float)*rvir
    vecs = np.identity(3,dtype = np.float)
    # Calculates eigensystem of intertia tensor
    inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))
    semmiaxes.append(np.copy(axes))
    eigenvecs.append(np.copy(vecs))
    # Process is repeated for R = 0.01*rvir
    axes = np.ones(3,dtype = np.float)*0.01*rvir
    vecs = np.identity(3,dtype = np.float)
    inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(vecs.ctypes.data), ct.c_void_p(axes.ctypes.data))
    semmiaxes.append(np.copy(axes))
    eigenvecs.append(np.copy(vecs))

# Reshapes so that each halo info is in the same row
semmiaxes = np.reshape(np.array(semmiaxes),(len(halonums),6))
eigenvecs = np.reshape(np.array(eigenvecs),(len(halonums),18))
np.savetxt("../Plots/Triaxiality/"+lvl+"/Semiaxes_rvir_1e2rvir.csv",semmiaxes,delimiter=',')
np.savetxt("../Plots/Triaxiality/"+lvl+"/Eigenvecs_rvir_1e2rvir.csv",eigenvecs,delimiter=',')

'''
############################################################
#    Level 3 Innerskirts DM vs MHD
############################################################
# Fonts 
# Fonts 
MEDIUM_SIZE = 30
SMALL_SIZE = 30
SSSMALL_SIZE = 16
plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        
# Simiulation specs
lvl3 = 'level3_MHD'
lvl3no = 'level3_DM'
lvl = 'level4_MHD'
lvlno = 'level4_DM'
# Obtains the axes
axes = np.loadtxt("../Plots/"+lvl+"/Semiaxes_rvir_1e2rvir.csv", delimiter = ',')#[(5,15,20,22,23,26),:]
axesno = np.loadtxt("../Plots/"+lvlno+"/Semiaxes_rvir_1e2rvir.csv", delimiter = ',')#[(5,15,20,22,23,26),:]
axes3 = np.loadtxt("../Plots/"+lvl3+"/Semiaxes_rvir_1e2rvir.csv", delimiter = ',')
axes3no = np.loadtxt("../Plots/"+lvl3no+"/Semiaxes_rvir_1e2rvir.csv", delimiter = ',')
# Axes in inner regions 0.01*rvir ~ 1-5 kpc
axes1 = axes[:,3:]
axes1_3 = axes3[:,3:]
# Axes in outer regions rvir ~ 100-500 kpc
axes2 = axes[:,:3]
axes2_3 = axes3[:,:3]
# Axes in inner regions 0.01*rvir ~ 1-5 kpc
axes1no = axesno[:,3:]
axes1no_3 = axes3no[:,3:]
# Axes in outer regions rvir ~ 100-500 kpc
axes2no = axesno[:,:3]
axes2no_3 = axes3no[:,:3]

print( "001rvir3 = ", (axes1_3[:,0]*axes1_3[:,1]*axes1_3[:,2])**(1./3.))
print( "001rvir4 = ", (axes1[:,0]*axes1[:,1]*axes1[:,2])**(1./3.))
print( "rvir3 = ", (axes2_3[:,0]*axes2_3[:,1]*axes2_3[:,2])**(1./3.))
print( "rvir4 = ", (axes2[:,0]*axes2[:,1]*axes2[:,2])**(1./3.))

# Plots axial ratios c/a Vs b/a for R = 0.01Rvir 
plt.plot(axes1_3[:,1]/axes1_3[:,0],axes1_3[:,2]/axes1_3[:,0], marker = 'o',
markersize =12, c = 'b', alpha = 0.8, linewidth = 0, label = r"$R_{MHD} = 0.01R_{200}$" )
plt.plot(axes1no_3[:,1]/axes1no_3[:,0],axes1no_3[:,2]/axes1no_3[:,0], marker = 's',
markersize =12, c = 'r', alpha = 0.8, linewidth = 0, label = r"$R_{DM} = 0.01R_{200}$" )


# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =12, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =12,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =12, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =12, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 3 Innerskirts DM vs MHD")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_Inner_lvl3.png",bbox_inches="tight")
plt.clf()







############################################################
#    Level 3 Outterskirts DM vs MHD
############################################################
# Fonts 
plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels


# Plots axial ratios c/a Vs b/a for R = Rvir 
plt.plot(axes2_3[:,1]/axes2_3[:,0],axes2_3[:,2]/axes2_3[:,0], marker = 'o',
c = 'blue',markersize =12, alpha = 0.8, linewidth = 0, label = r"$R_{MHD} = R_{200}$" )
plt.plot(axes2no_3[:,1]/axes2no_3[:,0],axes2no_3[:,2]/axes2no_3[:,0], marker = 's',
c = 'red',markersize =12, alpha = 0.8, linewidth = 0, label = r"$R_{DM} = R_{200}$" )


# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =12, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =12,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =12, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =12, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 3 Outterskirts DM vs MHD")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_Outter_lvl3.png",bbox_inches="tight")
plt.clf()






############################################################
#    Level 4 Innerskirts DM vs MHD
############################################################
# Fonts 

plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('figure', figsize=(7, 10))

plt.figure(figsize=(7, 10))
# Plots axial ratios c/a Vs b/a for R = 0.01Rvir 
plt.plot(axes1[:,1]/axes1[:,0],axes1[:,2]/axes1[:,0], marker = 's',
 c = 'r',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{MHD} = 0.01R_{200}$" )
plt.plot(axes1no[:,1]/axes1no[:,0],axes1no[:,2]/axes1no[:,0], marker = 'o',
c = 'b',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{DM} = 0.01R_{200}$" )

# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =15, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =15,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =15, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =15, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 4 Innerskirts DM vs MHD")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_Inner_lvl4.png",bbox_inches="tight")
plt.clf()




############################################################
#    Level 4 Outterskirts DM vs MHD
############################################################
# Fonts 

plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('figure', figsize=(7, 10))

plt.figure(figsize=(7, 10))
# Plots axial ratios c/a Vs b/a for R = Rvir 
plt.plot(axes2[:,1]/axes2[:,0],axes2[:,2]/axes2[:,0], marker = 's', 
 c = 'red',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{MHD} = R_{200}$" )
plt.plot(axes2no[:,1]/axes2no[:,0],axes2no[:,2]/axes2no[:,0], marker = 'o', 
c = 'blue',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{DM} = R_{200}$" )

# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =15, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =15,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =15, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =15, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 4 Innerskirts DM vs MHD")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_Outter_lvl4.png",bbox_inches="tight")
plt.clf()


############################################################
#    Level 4 DM Inner vs Outter
############################################################
# Fonts 

plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('figure', figsize=(7, 10))

plt.figure(figsize=(7, 10))
# Plots axial ratios c/a Vs b/a for R = Rvir and R= 0.01Rvir
plt.plot(axes1no[:,1]/axes1no[:,0],axes1no[:,2]/axes1no[:,0], marker = 's',
c = 'r',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{DM} = 0.01R_{200}$" )
plt.plot(axes2no[:,1]/axes2no[:,0],axes2no[:,2]/axes2no[:,0], marker = 'o', 
c = 'b',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{DM} = R_{200}$" )

# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =15, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =15,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =15, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =15, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 4 DM Inner vs Outter")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_DM_lvl4.png",bbox_inches="tight")
plt.clf()



############################################################
#    Level 4 MHD Inner vs Outter
############################################################
# Fonts 

plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('figure', figsize=(7, 10))

plt.figure(figsize=(7, 10))

# Plots axial ratios c/a Vs b/a for R = Rvir and R= 0.01Rvir
plt.plot(axes1[:,1]/axes1[:,0],axes1[:,2]/axes1[:,0], marker = 's', 
c = 'r',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{MHD} = 0.01R_{200}$" )
plt.plot(axes2[:,1]/axes2[:,0],axes2[:,2]/axes2[:,0], marker = 'o', 
c = 'b',markersize =12, alpha = 0.8, linewidth = 0, label = r"$(OUR) R_{MHD} = R_{200}$" )

# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =15, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =15,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =15, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =15, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 4 MHD Inner vs Outter")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_MHD_lvl4.png",bbox_inches="tight")
plt.clf()



############################################################
#    Level 3 DM Inner vs Outter
############################################################
# Fonts 
plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

# Plots axial ratios c/a Vs b/a for R = Rvir and R= 0.01Rvir
plt.plot(axes1no_3[:,1]/axes1no_3[:,0],axes1no_3[:,2]/axes1no_3[:,0], marker = 's',
c = 'r',markersize =12, alpha = 0.8, linewidth = 0, label = r"$R_{DM} = 0.01R_{200}$" )
plt.plot(axes2no_3[:,1]/axes2no_3[:,0],axes2no_3[:,2]/axes2no_3[:,0], marker = 'o',
c = 'b',markersize =12, alpha = 0.8, linewidth = 0, label = r"$R_{DM} = R_{200}$" )

# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =12, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =12,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =12, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =12, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 3 DM Inner vs Outter")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_DM_lvl3.png",bbox_inches="tight")
plt.clf()



############################################################
#    Level 3 MHD Inner vs Outter
############################################################
# Fonts 
plt.rc('font', size=SSSMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

# Plots axial ratios c/a Vs b/a for R = Rvir and R= 0.01Rvir
plt.plot(axes1_3[:,1]/axes1_3[:,0],axes1_3[:,2]/axes1_3[:,0], marker = 's',
c = 'r',markersize =12, alpha = 0.8, linewidth = 0, label = r"$R_{MHD} = 0.01R_{200}$" )
plt.plot(axes2_3[:,1]/axes2_3[:,0],axes2_3[:,2]/axes2_3[:,0], marker = 'o',
c = 'blue',markersize =12, alpha = 0.8, linewidth = 0, label = r"$R_{MHD} = R_{200}$" )

# Plots Observational references
plt.errorbar([1],[0.47], yerr = 0.14, label = "Loebman et al. @20kpc",marker = 'H',markersize =12, color = 'black')
plt.plot([1],[0.9], label = "Vera-Ciro et al. ~<10kpc",marker = '*',linewidth = 0, markersize =12,c = 'm')#<~ 10kpc
plt.plot([0.9],[0.8], label = "Vera-Ciro et al. >>30kpc",marker = '8',linewidth = 0, markersize =12, c = 'g')#>> 30kpc
plt.plot([0.99],[0.72], label = "Law & Majewski 2010 ",marker = 'P',linewidth = 0, markersize =12, c = 'y')# Must be outerskirts
plt.plot([0,1],[0,1])

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel("b/a")
plt.ylabel("c/a")
plt.title("Level 3 MHD Inner vs Outter")
plt.legend(loc = 0)
plt.savefig("../Plots/Triaxiality/"+"/Triaxiality_MHD_lvl3.png",bbox_inches="tight")
plt.clf()


