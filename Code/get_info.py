from arepo import *

def loadSim(lvl,halo,snapnumDM,snapnumMHD):
    
    # Path of the simulation
    pathDM = "/hits/universe/GigaGalaxy/"+lvl+"_DM"+"/"+halo+"/output" 
    pathMHD = "/hits/universe/GigaGalaxy/"+lvl+"_MHD"+"/"+halo+"/output"
    # Loads groups to select the principal halo
    subSnDM = Subfind(pathDM,snapnumDM)
    subSnMHD = Subfind(pathMHD,snapnumMHD)

    # Loads the simulation
    snDM = Snapshot(pathDM,snapnumDM, parttype=[1], combineFiles=True, verbose=False)
    snMHD = Snapshot(pathMHD,snapnumMHD, parttype=[1], combineFiles=True, verbose=False)
    
    # DM Mass
    DMassDM = snDM.MassTable[1]
    DMassMHD = snMHD.MassTable[1]
    
    # r the radius of the halo (R_mean200)*1000 to obtain it in kpc
    rvirDM = 1000.*subSnDM.group.Group_R_TopHat200[0]
    rvirMHD = 1000.*subSnMHD.group.Group_R_TopHat200[0]
    
    MvirDM = 1000.*subSnDM.group.Group_M_TopHat200[0]
    MvirMHD = 1000.*subSnMHD.group.Group_M_TopHat200[0]
    
    
    # Gets the length of the main structure particles /// 0 - Principal halo, 1 - DM
    subCutDM = subSnDM.subhalo.SubhaloLenType[0,1]
    subCutMHD = subSnMHD.subhalo.SubhaloLenType[0,1]
    
    lista = [subCutDM/1e+6,subCutMHD/1e+6, DMassDM/1e-5, DMassMHD/1e-5,rvirDM,rvirMHD,MvirDM/1e+4,MvirMHD/1e+4]
    stri = halo+"&"
    for el in lista[:-1]:
        stri += '{:.3f}'.format(float(el)) +"&"
    stri += '{:.3f}'.format(float(lista[-1]))+"\\\\"    
    
    print stri

#lvl = 'level4'
lvl = 'level3'

for i in [6,16,21,23,24,27]:#range(1,31):
    
    halo = 'halo_'+str(i)    
    loadSim(lvl,halo,127,63)
    
    
    
    
