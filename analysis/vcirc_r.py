from readData.DataLoader import DataLoader
import numpy as np
import matplotlib.pyplot as plt
import sys

from analysis import analyze
import shared_data

def calc_m(cat, orig_coords, part_types, gal_pos):

    mass_tot = dict()
    r2_tot = dict()
    for pt in part_types:
        coords = orig_coords[pt] - gal_pos
        mass = cat[f'PartType{pt}/Masses']*1e10/cat.h

        r2_tot[pt] = np.sum(np.square(coords), axis=1)
        mass_tot[pt] = mass

    return mass_tot, r2_tot

def main():

    path = "/home/j.rose/gr/"
    part_types = [4,0,1]
    keys = ["Coordinates", "Masses", "Velocities"] 

    G = 4.3e-6 #kpc/Msun * (km/s)^2

    path = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    #run = 'CDM_du10_s75_phys_no_formation_time/output-blue/'
    run = 'CDM_du10_s75_phys_correct_time/output-blue/'
    all_runs = [path+run]

    path = '/home/j.rose/Projects/SMUGGLE/RUNs/'
    run = 'CDM_du10/output-orange'
    all_runs += [path + run]

    snaps = [0, 75]
    names = ['Isolated', 'Zoom']

    fig, ax = shared_data.set_plot_params(ncols=3)
    for i,run in enumerate(all_runs):
        if 'isolated' in run:
            cat = DataLoader(run, part_types=part_types, snap_num=snaps[i], keys=keys)
            coords = dict()
            for pt in part_types:
                coords[pt] = cat[f'PartType{pt}/Coordinates'] / cat.h
        else:
            cat = DataLoader(run, part_types=part_types, snap_num=snaps[i], keys=keys+['SubhaloPos'], sub_idx=0)
            coords = dict()
            for pt in part_types:
                coords[pt] = cat[f'PartType{pt}/Coordinates'] / cat.h * cat.time

        if 'SubhaloPos' in cat:
            gal_pos = cat['SubhaloPos'] * cat.time / cat.h
        else:
            gal_pos = np.array([325]*3) / cat.h
        
        rads = [150,100,75,50,25]
        for rad in rads:
            box_cut = analyze.get_box_cut(coords[4], gal_pos, rad)
            gal_pos = np.average(coords[4][box_cut],axis=0,weights=cat['Masses'][box_cut])

        all_r = np.linspace(0.1, 1, 101)
        vcirc = {i:[] for i in [1,0,4,-1]}
        mass_type, r2_type = calc_m(cat, coords, part_types, gal_pos) 
        for r in all_r:
            for pt in part_types:
                 vcirc[pt].append(np.sqrt(G*np.sum(mass_type[pt][r2_type[pt] < r*r]) / r))
            vcirc[-1].append(np.sqrt(G*np.sum([np.sum(mass_type[pt][r2_type[pt] < r*r]) for pt in part_types]) /r))

        ax[i].plot(all_r, vcirc[-1], 'k-', label='total')
        ax[i].plot(all_r, vcirc[4], 'g:', label='stars')
        ax[i].plot(all_r, vcirc[0], 'b--', label='gas')
        ax[i].plot(all_r, vcirc[1], 'r-.', label='halo')

        p = ax[2].plot(all_r, vcirc[-1], '-', label=names[i])
        ax[2].plot(all_r, vcirc[4], ':', color=p[0].get_color())
        ax[2].plot(all_r, vcirc[0], '--', color=p[0].get_color())
        ax[2].plot(all_r, vcirc[1], '-.', color=p[0].get_color())

        ax[i].set_ylim((0,350))
        #ax[i//4,i%4].set_ylim((0,100))
        #ax[i//4,i%4].set_xlim((0,25))
        ax[i].set_xlabel("R [kpc]")
        ax[i].set_ylabel("$\mathrm{V}_{\mathrm{c}}$ [$\mathrm{kms}^{-1}$]")
        ax[i].set_title(names[i])
        

    ax[i].legend()

    ax[2].set_ylim((0,350))
    #ax[1,3].set_ylim((0,100))
    #ax[1,3].set_xlim((0,25))
    ax[2].set_xlabel("R [kpc]")
    ax[2].set_ylabel("$\mathrm{V}_{\mathrm{c}}$ [$\mathrm{kms}^{-1}$]")
    ax[2].legend()

    fig.savefig(f"plots/vcirc_r.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()
