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

    if type(fov) == type([]):
        xrange=[-fov[0]/2,fov[0]/2]
        yrange=[-fov[1]/2,fov[1]/2]
        maxfov = np.max(fov)
    else:
        xrange=[-fov/2,fov/2]
        yrange=[-fov/2,fov/2]
        maxfov = fov

    pt = part_types[0]

    coords = copy(cat[f'PartType{pt}/Coordinates'])
    vels = copy(cat[f'PartType{pt}/Velocities'])
    masses = copy(cat[f'PartType{pt}/Masses']) * 1e10

    data = dict()
    for key in map_keys:
        data[key] = copy(cat[f'PartType{pt}/{key}']) 
        if 'Mass' in key:
            data[key] *= 1e10
        if 'Density' in key and 'isolated' not in cat.path:
            data[key] /= cat.time**3
        if 'Pressure' in key and 'isolated' not in cat.path:
            data[key] /= cat.time**3 #check this
        if 'Softening' in key and 'isolated' not in cat.path:
            data[key] *= 1 #cat.time #check this
        print(key, np.min(data[key]), np.max(data[key]))

        print(np.sum(data[key] > 0))
    #return

    if 'Potential' in data:
        data['Potential'] *= -1

    coords -= gal_pos
    vels -= gal_vel

    if 'isolated' not in cat.path:
        coords *= cat.time
        vels *= np.sqrt(cat.time) 

    box_cut = analyze.get_box_cut(coords, np.zeros(3), maxfov)
    #box_cut = np.ones(len(coords), dtype=bool)

    coords = coords[box_cut]
    vels = vels[box_cut]
    masses = masses[box_cut]
    
    for key in map_keys:
        data[key] = data[key][box_cut]
        print(key, np.min(data[key]), np.max(data[key]))
        print(np.sum(data[key] > 0))

    coords, vels = analyze.get_rotate_data(coords, vels, masses, face_on=True, edge_on=False, r_min=1, r_max=7) 
    #phi, theta = analyze.get_rotate_data(coords, vels, masses, face_on=face_on, edge_on=edge_on, r_min=1, r_max=7, get_pt=True) 
    #print(phi, theta)
    #return
    #coords, vels = analyze.get_rotate_data(coords, vels, masses, phi=1.07, theta=2.15) 

    #cut = (coords[:,2] < 5) & (coords[:,2]>-5)
    cut = np.ones(coords.shape[0], dtype=bool)

    if 'SubfindHsml' in cat.data:
        hsml = cat['SubfindHsml']
    else:
        hsml = calc_hsml.get_particle_hsml(coords[cut,0], coords[cut,1], coords[cut,2])

    massmap,image = cmakepic.simple_makepic(coords[cut,0], coords[cut,1],
            weights = masses[cut], hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)

    maps = []
    for key in map_keys:
        print(key)

        if np.sum(data[key][cut] > 0) == 0:
            maps.append(np.ones((256,256)))
            continue

        if key == 'Velocities':
            map,image = cmakepic.simple_makepic(coords[cut,0], coords[cut,1],
                    weights = abs(vels[cut,0])*masses[cut], hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
        elif key == "GFM_Metals":
            for i in range(10):
                map,image = cmakepic.simple_makepic(coords[cut,0], coords[cut,1],
                        weights = data[key][:,i][cut]*masses[cut], hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
                maps.append(map / massmap)
        else:
            print()
            print(np.sum(data[key] > 0))
            weights = data[key][cut]*masses[cut]
            print(np.sum(weights > 0))
            print()
            map,image = cmakepic.simple_makepic(coords[cut,0], coords[cut,1],
                    weights = weights , hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
            print(np.sum(map > 0))

        if key != "GFM_Metals":
            maps.append(map / massmap)

        print(np.sum(maps[-1] > 0))

    return maps 

def main():

    pixels = 512 #can be made larger #does this actually do anything?
    fov = 20 #kpc - diameter (not radius)
    do_faces = True

    ipath = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    zpath = '/home/j.rose/Projects/SMUGGLE/RUNs/'
    keys = ["Coordinates", "Velocities", "Masses", "GFM_Metallicity", "StarFormationRate", "ElectronAbundance", 'InternalEnergy', "GFM_StellarFormationTime", "GFM_InitialMass", "VirialParameter", "Density", "GFM_Metals", 'Potential', 'Pressure', 'Softenings', 'GasRadCoolShutoffTime'] 

    zoom_runs = ['CDM_du10/output-orange/'] 
    all_runs = [zpath + run for run in zoom_runs]

    iso_runs = ['du10/output-orange/'] * 2
    all_runs += [ipath + run for run in iso_runs]

    iso_runs = ['du10_time/output-blue/'] * 5
    #all_runs += [ipath + run for run in iso_runs]

    ic_runs = ['../analysis/plots/ICs/'] * 2
    #all_runs += [ipath + run for run in ic_runs]

    snaps = [75, 0, 1, 2, 3, 4, 5, 6]
    names = ['Zoom', 'Isolated - S0', 'Isolated - S1', 'Isolated - S2', 'Isolated - S3', 'Isolated - S4', 'Isolated - S5']

    fig,ax = shared_data.set_plot_params(nrows=len(all_runs), ncols=17) #17 #12
    if len(ax.shape) == 1:
        ax = np.array([ax])

    for i,run in enumerate(all_runs):
        print(run)
        ax_idx = 0
        snap_num = snaps[i]

        if 'isolated' not in run:
            cat = DataLoader(run, part_types=[0,4], snap_num=snap_num, keys=keys+['GroupPos','GroupCM'], fof_idx=0) 
            gal_pos = cat['GroupPos']
            print(gal_pos)

        else:
            cat = DataLoader(run, part_types=[0,4], snap_num=snap_num, keys=keys)  
            gal_pos = np.average(cat['PartType4/Coordinates'],axis=0)
            gal_vel = np.average(cat['PartType4/Velocities'],axis=0)

        rads = [150,100,75,50,25]
        for rad in rads:
            box_cut = analyze.get_box_cut(cat['PartType4/Coordinates'], gal_pos, rad)
            gal_pos = np.average(cat['PartType4/Coordinates'][box_cut],axis=0,weights=cat['PartType4/Masses'][box_cut])
        gal_vel = np.average(cat['PartType4/Velocities'][box_cut],axis=0,weights=cat['PartType4/Masses'][box_cut])

        print(gal_pos)
        if 'isolated' in run:
            gal_pos = np.array([325,325,325])

        #PartType0 - Gas
        #map_keys = ["Velocities", "Masses", "GFM_Metallicity", "StarFormationRate", "ElectronAbundance", 'InternalEnergy'] #, "VirialParameter"]
        map_keys = ['Density', 'ElectronAbundance', 'GFM_Metallicity', 'GasRadCoolShutoffTime', 'InternalEnergy', 'Masses', 'Potential', 'Pressure', 'Softenings', 'StarFormationRate', 'Velocities', 'VirialParameter']
        #map_keys = ["GFM_Metals"]
        #map_keys = ['GasRadCoolShutoffTime']

        if len(map_keys) > 0:
            massmaps = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[0], filt=f'{names[i]}', gal_pos=gal_pos, gal_vel=gal_vel, map_keys=map_keys)
        else:
            massmaps = []

        vmin = [1e-6, 5e-3, 5e-3, 1e-14, 1e2, 2e4, 5e4, 1e-3, 2e-1, 1e-13, 5e0, 5e1]
        vmax = [1e0 , 1e0 , 7e-2, 1e-6 , 5e4, 5e4, 5e5, 1e2 , 2e0, 1e-3 , 2e2, 1e5]
        #vmin = [5e0, 1e4, 5e-3, 1e-8, 1e-2, 1e2, 5e1]
        #vmax = [3e2, 1e5, 1e-1, 1e-3, 1e0 , 1e4, 5e4]
        #vmin = [1e5] #volume
        #vmax = [1e11] #volume
        #vmin = [1e-5, 6e-6] + [1e-8]*8 #metals
        #vmax = [4e-5, 2e-5] + [1e-6]*8 #metals
        #titles = ['abs X Vel', 'Cell Mass', 'Metallicity', 'SFR', 'Electron Abundance', 'Internal Energy', "VirialParameter"]
        titles = ['Density', 'ElectronAbundance', 'Metallicity', 'GasRadCoolShutoffTime', 'InteralEnergy', 'Masses', 'Potential', 'Pressure', 'Softenings', 'StarFormationRate', 'Velocities', 'VirialParameter']
        #titles = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe', 'Other']
        for j,massmap in enumerate(massmaps):
            print("###", np.min(massmap), "###")
            im = ax[i,ax_idx].imshow(massmap, norm=LogNorm(vmin=vmin[j], vmax=vmax[j]))
            fig.colorbar(im, ax=ax[i,ax_idx])
            ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_xlabel('X [kpc]')
            ax[i,ax_idx].set_ylabel('Y [kpc]')

            ax[i,ax_idx].set_title(f"{names[i]} Gas - {titles[j]}")
            ax_idx += 1


        #PartType4 - stars
        map_keys = ["Velocities", "Masses", "GFM_Metallicity", "GFM_StellarFormationTime", "GFM_InitialMass"]
        #map_keys = ['GFM_StellarFormationTime']
        #map_keys = []
        if len(map_keys) > 0:
            massmaps = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[4], filt=f'{names[i]}', gal_pos=gal_pos, gal_vel=gal_vel, map_keys=map_keys)
        else:
            massmaps = []

        vmin = [2e1, 1e4, 1e-3, 2e0, 1e4]
        vmax = [2e2, 1e5, 2e-1, 5e0, 1e5]
        titles = ['abs X Vel', 'Particle Mass', 'Metallicity', 'Formation Time', 'Particle Initial Mass']
        #vmin = [2]
        #vmax = [5]
        #titles = ['Formation Time']
        for j,massmap in enumerate(massmaps):
            print("###", np.min(massmap), "###")
            #if i in [1,3]:
            #    massmap = 13.7 - units.age_from_a(massmap, H0=69.09, Om0=0.301712)
            im = ax[i,ax_idx].imshow(massmap, norm=LogNorm(vmin=vmin[j], vmax=vmax[j]))
            fig.colorbar(im, ax=ax[i,ax_idx])
            ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_xlabel('X [kpc]')
            ax[i,ax_idx].set_ylabel('Y [kpc]')

            ax[i,ax_idx].set_title(f"{names[i]} Stars - {titles[j]}")
            ax_idx += 1

    fig.savefig(f"plots/light_map_{snap_num}.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()

