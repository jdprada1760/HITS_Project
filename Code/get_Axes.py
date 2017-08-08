import matplotlib
matplotlib.use('Agg')

# Python
from arepo import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import yt
from yt.units import parsec, Msun
from Py_Libs.shape import *

# C imports
import ctypes as ct
inertia = ct.cdll.LoadLibrary('./C_Libs/inertia.so')

# Simiulation specs
#lvl = 'level5'
lvl = 'level4_MHD'
#lvl = 'level5_Durham' 
halonums = [20]
#halonums = range(1,11)
#halonums = ['halo16','halo16_MHD','halo24','halo24_MHD','halo28', 'halo6_MHD','halo9','halo9_MHD']

# Loads simulation
snapnum = 127
pos,rvir = loadSim(lvl,'halo_1',snapnum)
pos = np.array(pos,dtype = np.float)
eigenvals = rvir*np.ones(3,dtype = np.float)
eigenvecs = np.identity(3,dtype = np.float)
inertia.get_shape(ct.c_void_p(pos.ctypes.data),ct.c_int(len(pos)), ct.c_void_p(eigenvecs.ctypes.data), ct.c_void_p(eigenvals.ctypes.data))
print(eigenvals)
print(eigenvecs)
