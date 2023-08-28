import matplotlib
matplotlib.use("agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys, os, time
from copy import copy

from readData.DataLoader import DataLoader 
from analysis import analyze

import units.springel_units as units  
from scipy.stats import binned_statistic

import visualization.contour_makepic as cmakepic
import util.calc_hsml as calc_hsml
import shared_data

from mpi4py import MPI

def mag_to_lum(V, mag_sun):

    lums = np.zeros(V.shape[0])
    for i in range(V.shape[0]):
        mag = V[i]
        lums[i] = np.power(10, (mag - mag_sun)/2.5*-1)
    return lums

def get_massmap(cat, pixels, fov=50, face_on=False, edge_on=False, part_types=[4], filt="V", is_light=True, gal_pos=None, gal_vel=None):

    if type(fov) == type([]):
        xrange=[-fov[0]/2,fov[0]/2]
        yrange=[-fov[1]/2,fov[1]/2]
        maxfov = np.max(fov)
    else:
        xrange=[-fov/2,fov/2]
        yrange=[-fov/2,fov/2]
        maxfov = fov

    coords = cat['Coordinates'] #*cat.time
    vels = cat['Velocities']
    
    masstab = cat.masstable[part_types[0]]       
    if masstab == 0:
        masses = cat['Masses']
    else:
        masses = np.array([masstab]*len(vels))    

    if np.max(np.max(coords, axis=0) - np.min(coords, axis=0) > 600):
        gal_pos += 325
        coords += 325

        gal_pos -= 650
        coords -= 650

    coords -= gal_pos
    vels -= gal_vel

    box_cut = analyze.get_box_cut(coords, np.zeros(3), maxfov)
    box_cut = np.ones(len(coords), dtype=bool)

    coords = coords[box_cut]
    vels = vels[box_cut]
    masses = masses[box_cut]

    if 'SubfindHsml' in cat.data:
        hsml = cat['SubfindHsml']
    else:
        hsml = calc_hsml.get_particle_hsml(coords[:,0], coords[:,1], coords[:,2])

    face_coords, face_vels = analyze.get_rotate_data(copy(coords), copy(vels), masses, phi=0, theta=0)  

    face_massmap,image = cmakepic.simple_makepic(face_coords[:,0], face_coords[:,1],
            weights=masses, hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)

    return face_massmap

def main():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    pixels = 256
    fov = 300
    snap_nums = 500
    num_snaps = 500
    skip = 1

    path = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    keys = ["Coordinates", "Velocities", "Masses"] 

    run = 'du10_large/output-blue/'
    name = 'du10-large-dm'

    if rank == 0:
        if not os.path.exists(f"plots/{name}/"):
            os.mkdir(f"plots/{name}")
    else:
        time.sleep(.2)

    numbers = np.arange(snap_nums, snap_nums+num_snaps, skip)
    indexes = np.where(np.arange(numbers.shape[0])%size==rank)[0]
    realizations = numbers[indexes]

    print(rank, realizations)

    for j,snap_num in enumerate(realizations):

        fig,ax = shared_data.set_plot_params()
    
        cat = DataLoader(path+run, part_types=[1], snap_num=snap_num, keys=keys) 

        gal_pos = np.array([325]*3)
        rads = [325,150,100,75,50] #,25]
        for rad in rads:
            box_cut = analyze.get_box_cut(cat['Coordinates'], gal_pos, rad)
            gal_pos = np.average(cat['Coordinates'][box_cut],axis=0,weights=cat['Masses'][box_cut])
        gal_vel = np.average(cat['Velocities'][box_cut],axis=0,weights=cat['Masses'][box_cut])
        print(snap_num, gal_pos)
        print(gal_vel)

        #gal_pos = np.array([325,325,325])
        gal_pos_orig = copy(gal_pos)

        face_massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[1], filt='', gal_pos=gal_pos, gal_vel=gal_vel)

        im = ax.imshow(face_massmap, norm=LogNorm(vmin=5e-4,vmax=2e-1))
        #fig.colorbar(im, ax=ax)
        ax.set_xticks([pixels/5*j for j in range(6)])
        ax.set_yticks([pixels/5*j for j in range(6)])
        ax.set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
        ax.set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
        ax.set_xlabel("X [kpc]")
        ax.set_ylabel("Y [kpc]")

        ax.set_title(name)

        fig.savefig(f"plots/{name}/light_map_{snap_num:04d}.png", dpi=100, bbox_inches='tight')
        print(f"plots/{name}/light_map_{snap_num:04d}.png")
        plt.close(fig)

    return

if __name__=="__main__":
    main()
