import matplotlib
matplotlib.use('Agg')

from arepo import *
import numpy as np
import matplotlib.pyplot as plt
import yt
from yt.units import parsec, Msun

# Specification of the simulation 
path = "/hits/universe/GigaGalaxy/level4_MHD_new/halo_L2/output/snapdir_127" 
#snap = "/snapshot_063.hdf5"
#subSnap = "/fof_subhalo_tab_063.hdf5"

ds = yt.load(path)


##############################################################################################
# VOLUME RENDERING
##############################################################################################
field = 'gas','density'
im, sc = yt.volume_render(ds, field=field, fname="v0.png", sigma_clip=6.0)

sc.camera.set_width(ds.arr(0.1,'code_length'))
tf = sc.get_source(0).transfer_function
tf.clear()
tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
        alpha=np.logspace(-3,0,4), colormap = 'RdBu_r')
sc.render()
sc.save("v1.png", sigma_clip=6.0)
