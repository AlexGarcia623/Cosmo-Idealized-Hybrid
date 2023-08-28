import matplotlib
matplotlib.use("agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys, os
from copy import copy

from readData.DataLoader import DataLoader 
from analysis import analyze

import visualization.contour_makepic as cmakepic
import util.calc_hsml as calc_hsml
import shared_data

def mag_to_lum(V, mag_sun):

    lums = np.zeros(V.shape[0])
    for i in range(V.shape[0]):
        mag = V[i]
        lums[i] = np.power(10, (mag - mag_sun)/2.5*-1)
    return lums

def get_massmap(cat, pixels, fov=50, face_on=False, edge_on=False, part_types=[4], filt="V", name=''):

    if type(fov) == type([]):
        xrange=[-fov[0],fov[0]]
        yrange=[-fov[1],fov[1]]
        maxfov = np.max(fov)
    else:
        xrange=[-fov,fov]
        yrange=[-fov,fov]
        maxfov = fov

    coords = copy(cat['Coordinates'])/cat.h # *cat.time
    vels = copy(cat['Velocities'])
    masses = cat['Masses']*1e10/cat.h

    #gal_pos = cat['SubhaloPos']/cat.h # *cat.time
    #gal_vel = cat['SubhaloVel']

    gal_pos = cat['GroupPos']/cat.h # *cat.time
    gal_vel = cat['GroupVel']

    #info from Willmer et al. 2018 (K is Kshort - not sure if this is right)
    filt_idx = {'U':0, 'B':1, 'V':2, "K":3, "g":4, "r":5, "i":6, "z":7}
    solar_correction = {'U':6.33, "B":5.31, "V":4.80, "K":6.39, "g":5.11, "r":4.61, "i":4.52, "z":4.50}

    if 4 in part_types:
        phots = cat['GFM_StellarPhotometrics']
        lums = mag_to_lum(phots[:, filt_idx[filt]], solar_correction[filt]) #convert to luminosity 4.80=V (AB abs)
        weights = lums
    else:
        weights = masses

    box_cut = analyze.get_box_cut(coords, gal_pos, maxfov)

    coords -= gal_pos
    vels -= gal_vel

    coords, vels = analyze.get_rotate_data(coords, vels, masses, face_on=face_on, edge_on=edge_on) 

    if 4 in part_types:
        zcut = (coords[:,2] < 10) & (coords[:,2]>-10)
        cut = zcut & box_cut
    else:
        cut = box_cut

    if 'PartType1/SubfindHsml' in cat.data:
        hsml = cat['SubfindHsml']
    elif 'PartType4/SubfindHsml' in cat.data:
        hsml = cat['SubfindHsml']
    #elif f"{name}.npy" in os.listdir("np_arrays/"):
    #    hsml = np.load(f"np_arrays/{name}.npy")
    else:
        hsml = calc_hsml.get_particle_hsml(coords[:,0], coords[:,1], coords[:,2], DesNgb=128)
        #np.save(f"np_arrays/{name}", hsml)

    massmap,image = cmakepic.simple_makepic(coords[cut,0], coords[cut,1],
            weights=weights[cut], hsml=hsml[cut],xrange=xrange,yrange=yrange, pixels=pixels)

    return massmap


def main():

    path = "/home/j.rose/Projects/SMUGGLE/isolated/RUNs/"
    part_types = [1]
    keys = ["Coordinates", "Velocities", "Masses", "SubhaloPos", "GFM_StellarPhotometrics", "SubhaloVel", 'SubfindHsml', 'GroupPos', 'GroupVel'] 
    all_runs = ['1keV','30keV']
    snap_num = 51

    pixels = 512 

    sub_idx = -1
    fof_idx = 0

    filters = ['i','r','g']

    face_brightness = [0,0,0]
    edge_brightness = [0,0,0]

    plt.rc('font',**{'family':'STIXGeneral'})
    plt.rc('text', usetex=True)

    plt.rc('font', size=12)          # controls default text sizes
    plt.rc('axes', titlesize=16)     # fontsize of the axes title
    plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
    plt.rc('legend', fontsize=14)

    for i,run in enumerate(all_runs):
        print(run)

        cat = DataLoader(path+run, part_types=part_types, snap_num=snap_num, keys=keys, sub_idx=sub_idx, fof_idx=fof_idx) 

        if 0 in part_types:
            for side in ['face','edge']:

                if side=='face':
                    fig, ax = plt.subplots(figsize=(21,21))
                else:
                    fig, ax = plt.subplots(figsize=(21,7))

                if side == 'face':
                    fov = 15
                    massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=part_types, name=names[i])
                else:
                    fov = [5,15]
                    massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=True, part_types=part_types, name=names[i])


                ma = 9
                mi = 6
                massmap = np.log(massmap / mi) / np.log(ma/ mi)

                im = ax.imshow(massmap) 
                fig.colorbar(im, ax=ax)
                ax.set_xticks([])
                ax.set_yticks([])

                if side=='face':
                    fig.savefig(f"plots/gas_map_{names[i]}_face_{str(snap_num).zfill(3)}.pdf", bbox_inches='tight') #, dpi=200)
                else:
                    fig.savefig(f"plots/gas_map_{names[i]}_edge_{str(snap_num).zfill(3)}.pdf", bbox_inches='tight') #, dpi=200)

                plt.close(fig)

        if 1 in part_types or 2 in part_types:
            fig, ax = plt.subplots(figsize=(5,5))
            fov = 100
            massmap = get_massmap(cat, pixels=pixels, fov=fov, part_types=part_types, name=run) 
            massmap[np.isnan(massmap)] = 1e5
            massmap[massmap<=1e5] = 1e5

            im = ax.imshow(massmap, norm=LogNorm(vmin=1e5, vmax=1e8), cmap='gnuplot2')
            ax.set_xticks([])
            ax.set_yticks([])           

            fig.savefig(f"plots/mass_map_{run}.pdf", bbox_inches='tight')
            plt.close(fig)

        if 4 in part_types:
            for side in ['face','edge']:

                if side=='face':
                    fig, ax = plt.subplots(figsize=(21,21))
                else:
                    fig, ax = plt.subplots(figsize=(21,7))

                filter_maps = {filt:None for filt in filters}
                for filt in filters:
                    if side == 'face':
                        fov = 15
                        massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=part_types, filt=filt, name=names[i])
                    else:
                        fov = [5,15]
                        massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=True, part_types=part_types, filt=filt, name=names[i])
                    filter_maps[filt] = massmap 

                color_map = np.array([filter_maps[filt] for filt in filters])
                color_map = np.power(color_map, cat.time * .23) #.24) # + cat.redshift/5) #.316)
                shape = color_map.shape[1:]

                print(np.min(color_map), np.max(color_map), np.average(color_map), np.std(color_map))

                if side == 'face':
                    #mi = 32 #z0
                    #ma = 150 #z0
                    #mi = 20 #z0
                    #ma = 150 #z0
                    mi = 10
                    ma = 40
                    #if 'SMU' in names[i]:
                    #    mi = 3 #z1
                    #    ma = 10 #z1
                    #else:
                    #    mi = 5 #z1
                    #    ma = 13 #z1
                else:
                    #mi = 40 #z0 913
                    #ma = 180 #z0 913
                    #mi = 30 #z0
                    #ma = 180 #z0
                    mi = 10
                    ma = 40
                    #if 'SMU' in names[i]:
                    #    mi = 3 #z1
                    #    ma = 13 #z1
                    #else:
                    #    mi = 5 #z1
                    #    ma = 16 #z1
                #mi = np.average(color_map) #- np.std(color_map)
                #color_map[color_map < mi] = mi
                #color_map[color_map > ma] = ma
                for j in range(len(color_map)):
                    color_map[j] = np.log(color_map[j] / mi) / np.log(ma/ mi)

                #if side == 'face' and i == 0:
                #    face_brightness = [np.max(color_map[i]) for i in range(len(filters))]
                #elif side == 'edge' and i == 0:
                #    edge_brightness = [np.max(color_map[i]) for i in range(len(filters))]

                #if side == 'face':
                #    color_map[:3] /= np.array([np.ones(shape)*face_brightness[0],np.ones(shape)*face_brightness[1],np.ones(shape)*face_brightness[2]])
                #if side == 'edge':
                #    color_map[:3] /= np.array([np.ones(shape)*edge_brightness[0],np.ones(shape)*edge_brightness[1],np.ones(shape)*edge_brightness[2]])

                color_map = np.transpose(color_map, (1,2,0))           

                ax.imshow(color_map) 
                ax.set_xticks([])
                ax.set_yticks([])

                if side=='face':
                    fig.savefig(f"plots/light_map_{names[i]}_face_{str(snap_num).zfill(3)}.pdf", bbox_inches='tight') #, dpi=200)
                else:
                    fig.savefig(f"plots/light_map_{names[i]}_edge_{str(snap_num).zfill(3)}.pdf", bbox_inches='tight') #, dpi=200)
                continue

                if side=='face':
                    fig.savefig(f"plots/movies/{names[i]}/light_map_{names[i]}_face_{str(snap_num).zfill(3)}.pdf", bbox_inches='tight') #, dpi=200)
                else:
                    fig.savefig(f"plots/movies/{names[i]}/light_map_{names[i]}_edge_{str(snap_num).zfill(3)}.pdf", bbox_inches='tight') #, dpi=200)

                plt.close(fig)

    return

if __name__=="__main__":
    main()

