import matplotlib
matplotlib.use("agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys, os
from copy import copy

from readData.DataLoader import DataLoader 
from analysis import analyze
import units.springel_units as units  

import visualization.contour_makepic as cmakepic
import util.calc_hsml as calc_hsml
import shared_data

def mag_to_lum(V, mag_sun):

    lums = np.zeros(V.shape[0])
    for i in range(V.shape[0]):
        mag = V[i]
        lums[i] = np.power(10, (mag - mag_sun)/2.5*-1)
    return lums

def get_massmap(cat, pixels, fov=50, face_on=False, edge_on=False, part_types=[4], filt="V", is_light=True, gal_pos=None, gal_vel=None, map_keys=[]):


    pt = part_types[0]

    coords = copy(cat[f'PartType{pt}/Coordinates']) #*cat.time
    vels = copy(cat[f'PartType{pt}/Velocities'])
    masses = copy(cat[f'PartType{pt}/Masses']) * 1e10

    data = dict()
    for key in map_keys:
        data[key] = copy(cat[f'PartType{pt}/{key}']) 
        if 'Mass' in key:
            data[key] *= 1e10

    #data['Potential'] *= -1

    coords -= gal_pos
    vels -= gal_vel

    if 'isolated' not in cat.path:
        coords *= cat.time
        vels *= np.sqrt(cat.time) 

    box_cut = analyze.get_box_cut(coords, np.zeros(3), fov)
    box_cut = np.ones(len(coords), dtype=bool)

    coords = coords[box_cut]
    vels = vels[box_cut]
    masses = masses[box_cut]
   
    zcut = masses > 0
    print(np.sum(~zcut))

    for key in map_keys:
        data[key] = data[key][box_cut & zcut]

    maps = []
    for key in map_keys:
        mi = np.min(data[key])
        ma = np.max(data[key])
        avg = np.average(data[key])
        std = np.std(data[key])
        maps.append((key, mi, ma, avg, std))

    return maps 

def main():

    pixels = 512 #can be made larger #does this actually do anything?
    fov = 20 #kpc - diameter (not radius)
    do_faces = True

    ipath = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    zpath = '/home/j.rose/Projects/SMUGGLE/RUNs/'
    keys = ["Coordinates", "Velocities", "Masses", "GFM_Metallicity", "StarFormationRate", "ElectronAbundance", 'InternalEnergy', "GFM_StellarFormationTime", "GFM_InitialMass", "VirialParameter", "Density", "GFM_Metals"] 

    #iso_runs = ['new_vel/output-blue/']
    #iso_runs = ['new_arepo/output-blue/']
    iso_runs = ['new_order/output-blue/']
    all_runs = [ipath + run for run in iso_runs]

    zoom_runs = ['CDM_du10/output-orange/']
    all_runs += [zpath + run for run in zoom_runs]

    ic_runs = ['../analysis/plots/ICs/',] #'../analysis/plots/ICs/', '../analysis/plots/ICs/']
    all_runs += [ipath + run for run in ic_runs]

    snaps = [0, 75, 3] #, 1, 2]
    names = ['Iso ', 'Zoom', ' ICs ']

    gas_values = []
    star_values = []
    for i,run in enumerate(all_runs):
        print(run)
        ax_idx = 0
        snap_num = snaps[i]

        if 'isolated' not in run:
            cat = DataLoader(run, part_types=[0,4], snap_num=snap_num, keys=keys+['GroupPos','GroupCM'], fof_idx=0) 
            gal_pos = cat['GroupPos']

        else:
            cat = DataLoader(run, part_types=[0,4], snap_num=snap_num, keys=keys)  
            gal_pos = np.array([325,325,325])
            gal_vel = np.average(cat['PartType4/Velocities'],axis=0) 

        #PartType0 - Gas
        map_keys = ["Velocities", "Masses", "GFM_Metallicity", "StarFormationRate", "ElectronAbundance", 'InternalEnergy'] #, "VirialParameter"]
        #map_keys = ["GFM_Metals"]
        #map_keys = []

        if len(map_keys) > 0:
            massmaps = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[0], filt=f'{names[i]}', gal_pos=gal_pos, gal_vel=gal_vel, map_keys=map_keys)

        gas_values.append(massmaps)

        titles = ['abs X Vel', 'Cell Mass', 'Metallicity', 'SFR', 'Electron Abundance', 'Internal Energy', "VirialParameter"]
        #titles = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe', 'Other']


        #PartType4 - stars
        map_keys = ["Velocities", "Masses", "GFM_Metallicity", "GFM_StellarFormationTime", "GFM_InitialMass"]
        #map_keys = ['GFM_StellarFormationTime']
        #map_keys = []
        massmaps = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[4], filt=f'{names[i]}', gal_pos=gal_pos, gal_vel=gal_vel, map_keys=map_keys)

        star_values.append(massmaps)

        titles = ['abs X Vel', 'Particle Mass', 'Metallicity', 'Formation Time', 'Particle Initial Mass']


    print("Gas")
    for i in range(len(gas_values[0])):
        print(gas_values[0][i][0])
        for j in range(len(names)):
            key, mi, ma, avg, std = gas_values[j][i]
            print(f"{names[j]} - Min {mi:.2e} Max {ma:.2e} Avg {avg:.2e} Std {std:.2e}")
        print()

    print("Stars")
    for i in range(len(star_values[0])):
        print(star_values[0][i][0])
        for j in range(len(names)):
            key, mi, ma, avg, std = star_values[j][i]
            print(f"{names[j]} - Min {mi:.2e} Max {ma:.2e} Avg {avg:.2e} Std {std:.2e}")
        print()

    return

if __name__=="__main__":
    main()

