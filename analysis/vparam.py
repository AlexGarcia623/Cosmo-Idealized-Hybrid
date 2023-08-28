import matplotlib
matplotlib.use("agg")

from readData.DataLoader import DataLoader
import numpy as np
import matplotlib.pyplot as plt
import sys
import shared_data

import analysis.analyze as analyze
import os
import util.calc_hsml as calc_hsml
from scipy.stats import binned_statistic

def get_vparam(cat, r, gal_pos):

    #potential = cat['VirialParameter']
    potential = cat['Masses'] / cat['Density']
    coords = (cat['Coordinates'] - gal_pos) / cat.h

    if cat.time < 1:
        coords *= cat.time
        potential *= cat.time**3

    coords_r = np.sqrt(np.sum(np.square(coords[:,:2]),axis=1))
    potential_r = np.zeros(len(r))
    for i,rad in enumerate(r):
        rcut = (coords_r > rad) & (coords_r < rad+.2)
        potential_r[i] = np.median(potential[rcut])

    return potential_r

def main():

    path = "/home/j.rose/gr/"
    snap_num = 127
    part_types = [0]
    keys = ["Masses", "Velocities", "Coordinates", 'Potential', 'VirialParameter', 'Density'] 

    soft=.305/2
    ipath = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    zpath = '/home/j.rose/Projects/SMUGGLE/RUNs/'

    iso_runs = ['new_metals/output-blue/']
    all_runs = [ipath + run for run in iso_runs]

    zoom_runs = ['CDM_du10/output-orange/']
    all_runs += [zpath + run for run in zoom_runs]
    snaps = [0,75]
    fig, ax = shared_data.set_plot_params()

    alphas = [1]*10
    names = ['Isolated', 'Zoom']

    for i,run in enumerate(all_runs):

        snap_num = snaps[i]

        if 'isolated' not in run:
            cat = DataLoader(run, part_types=part_types, snap_num=snap_num, keys=keys+['GroupPos','GroupCM'], fof_idx=0) 
            gal_pos = cat['GroupPos']
        else:
            cat = DataLoader(run, part_types=part_types, snap_num=snap_num, keys=keys)  

        if cat.time > 1:
            gal_pos = np.array([325]*3)

        plot_r = np.logspace(-1,2,100)
        vparam = get_vparam(cat, plot_r, gal_pos)

        ax.plot(plot_r, vparam, label=names[i])

    ax.legend()
    ax.set_xlabel("Radius [kpc]")
    #ax.set_ylabel("Virial Parameter [?]")
    ax.set_ylabel("Cell Volume [kpc$^3$]")
    ax.set_xscale("log")
    ax.set_yscale("log")

    fig.savefig(f"plots/volume_profile.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()
