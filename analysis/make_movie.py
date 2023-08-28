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

        #cutg = gal_pos > 650
        #cutc = coords > 650
        #gal_pos[cutg] -= 650
        #coords[cutc] -= 650
        gal_pos -= 650
        coords -= 650

    coords -= gal_pos
    vels -= gal_vel

    box_cut = analyze.get_box_cut(coords, np.zeros(3), maxfov)
    box_cut = np.ones(len(coords), dtype=bool)

    coords = coords[box_cut]
    vels = vels[box_cut]
    masses = masses[box_cut]

    #edge_on = False
    #face_on = True
    #phi, theta = analyze.get_rotate_data(coords, vels, masses, face_on=face_on, edge_on=edge_on, r_min=1, r_max=7, get_pt=True) 
    #print(phi, theta)
    #return

    if 'SubfindHsml' in cat.data:
        hsml = cat['SubfindHsml']
    else:
        hsml = calc_hsml.get_particle_hsml(coords[:,0], coords[:,1], coords[:,2])

    face_coords, face_vels = analyze.get_rotate_data(copy(coords), copy(vels), masses, face_on=True)
    #face_coords, face_vels = analyze.get_rotate_data(copy(coords), copy(vels), masses, phi=-2.542, theta=1.367)  

    face_massmap,image = cmakepic.simple_makepic(face_coords[:,0], face_coords[:,1],
            weights=masses, hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
 

    edge_coords, edge_vels = analyze.get_rotate_data(coords, vels, masses, edge_on=True) 
    #edge_coords, edge_vels = analyze.get_rotate_data(coords, vels, masses, phi=-2.542, theta=2.938)

    edge_massmap,image = cmakepic.simple_makepic(edge_coords[:,0], edge_coords[:,1],
            weights=masses, hsml=hsml,xrange=[-3,3],yrange=yrange, pixels=pixels)

    return face_massmap, edge_massmap


def main():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    pixels = 512 
    fov = 15 #z2 25
    snap_nums = 0
    num_snaps = 24
    skip = 1

    path = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    keys = ["Coordinates", "Velocities", "Masses"] 


    run = 'du10_refine_highres/output-blue/'
    name = 'du10-refine-highres'

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

        fig,ax = shared_data.set_plot_params(ncols=3)
    
        cat = DataLoader(path+run, part_types=[4], snap_num=snap_num, keys=keys) 

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

        cat = DataLoader(path+run, part_types=[0], snap_num=snap_num, keys=keys) 

        face_massmap, edge_massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[0], filt='', gal_pos=gal_pos, gal_vel=gal_vel)

        im = ax[0].imshow(face_massmap, norm=LogNorm(vmin=5e-4,vmax=6e-1))
        #fig.colorbar(im, ax=ax[0])
        ax[0].set_xticks([pixels/5*j for j in range(6)])
        ax[0].set_yticks([pixels/5*j for j in range(6)])
        ax[0].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
        ax[0].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
        ax[0].set_xlabel("X [kpc]")
        ax[0].set_ylabel("Y [kpc]")

        ax[0].set_title(name)

        ax[1].imshow(edge_massmap, norm=LogNorm(vmin=5e-4, vmax=6e-1))
        ax[1].set_xticks([edge_massmap.shape[1]/5*j for j in range(6)])
        ax[1].set_yticks([edge_massmap.shape[0]/6*j for j in range(7)])
        ax[1].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
        ax[1].set_yticklabels([f'{j:.0f}' for j in range(7)])
        ax[1].set_xlabel("X [kpc]")
        ax[1].set_ylabel("Y [kpc]")

        cat = DataLoader(path+run, snap_num, part_types=4, keys=["Coordinates", "GFM_StellarFormationTime", "GFM_InitialMass", 'BirthPos']) 

        coords = cat['Coordinates']
        gal_pos = gal_pos_orig
        gal_r = np.sqrt(np.sum(np.square(gal_pos-np.array([325,325,325]))))

        if gal_r > 300:
            gal_pos += 325
            coords += 325

            gal_pos %= 325
            coords %= 325

        box_cut = analyze.get_box_cut(coords, gal_pos, fov)

        ages = cat['GFM_StellarFormationTime'][box_cut]
        masses = cat['GFM_InitialMass'][box_cut] * 1e10/cat.h
        cut = ages > 0 #wind particles
        if np.sum(cut) == 0:
            print("Fail")
            continue
        ages = ages[cut]
        masses = masses[cut]
        sfr_bins = np.linspace(5.6, 5.75, 200) #du10-short
        #sfr_bins = np.linspace(5.0, 7.5, 1000) #du10
        #sfr_bins = np.linspace(2.5, 6.5, 1000) #du10_z2
        #sfr_bins = np.linspace(9.0, 11.5, 1000) #du10_z.4
        mass_binned, bin_edges, bin_idx = binned_statistic(ages, masses, 'sum', bins=sfr_bins)
        bin_half_width = (bin_edges[1]-bin_edges[0])/2
        bins = bin_edges[:-1] + bin_half_width
        sfr = mass_binned / (bin_half_width*2*1e9)
        ax[2].plot(bins, sfr, label='Isolated')

        ax[2].set_xlabel("Look back time (Gyr)")
        ax[2].set_ylabel("SFR ($M_\odot$/yr)")
        ax[2].set_ylim((0, 25))

        #ages = np.load("plots/sfr.npy")
        #masses = np.load('plots/masses.npy')
        #mass_binned, bin_edges, bin_idx = binned_statistic(ages, masses, 'sum', bins=sfr_bins)
        #sfr = mass_binned / (bin_half_width*2*1e9)
        #ax[2].plot(bins, sfr, label='Zoom')
        #ax[2].legend(loc='upper right')


        fig.savefig(f"plots/{name}/light_map_{snap_num:04d}.png", dpi=100, bbox_inches='tight')
        print(f"plots/{name}/light_map_{snap_num:04d}.png")
        plt.close(fig)

    return

if __name__=="__main__":
    main()
