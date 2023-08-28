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

def get_param(cat, r, gal_pos, plot_param):

    param = cat[plot_param] #'VirialParameter'
    print(param.shape)

    if 'isolated' not in cat.path:
        coords = (cat['Coordinates'] - gal_pos) / cat.h * cat.time 
    else:
        coords = (cat['Coordinates'] - gal_pos) / cat.h  

    coords_r = np.sqrt(np.sum(np.square(coords[:,:2]),axis=1))
    param_r = np.zeros(len(r))
    for i,rad in enumerate(r):
        rcut = (coords_r > rad) & (coords_r < rad+.2)
        param_r[i] = np.median(param[rcut])
        #param_r[i] = np.sum(rcut)

    return param_r

def main():

    plot_param = 'Masses'

    ipath = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    zpath = '/home/j.rose/Projects/SMUGGLE/RUNs/'

    iso_runs = ['du10_refine_highres/output-blue/'] * 4
    all_runs = [ipath + run for run in iso_runs]

    zoom_runs = ['CDM_du10/output-orange/'] * 3
    all_runs += [zpath + run for run in zoom_runs]

    ic_runs = ['../analysis/plots/ICs/'] * 1
    all_runs += [ipath + run for run in ic_runs]

    snaps = [0,1,2,5,73,75,77,2,3,4,5,75,1]
    fig, ax = shared_data.set_plot_params()

    alphas = [1]*10
    names = [f'Isolated - S{snaps[i]}' for i in range(len(iso_runs))]
    names += ['Zoom 73', 'zoom 75', 'zoom 77','du10-no-taper']

    keys = ["Masses", "Velocities", "Coordinates"] #, plot_param] 
    part_types = [0,4]

    for i,run in enumerate(all_runs):

        snap_num = snaps[i]

        if 'isolated' not in run:
            cat = DataLoader(run, part_types=part_types, snap_num=snap_num, keys=keys+['GroupPos','GroupCM'], fof_idx=0) 
            gal_pos = cat['GroupPos']

        else:
            cat = DataLoader(run, part_types=part_types, snap_num=snap_num, keys=keys)  
            gal_pos = np.array([325]*3)

        if plot_param not in cat:
            continue

        if True:
            rads = [150,100,75,50,25]
            for rad in rads:
                box_cut = analyze.get_box_cut(cat['PartType4/Coordinates'], gal_pos, rad)
                gal_pos = np.average(cat['PartType4/Coordinates'][box_cut],axis=0,weights=cat['PartType4/Masses'][box_cut])
            gal_vel = np.average(cat['PartType4/Velocities'][box_cut],axis=0,weights=cat['PartType4/Masses'][box_cut])
        else:
            if 'isolated' in run:
                gal_pos = np.array([325]*3)

        #do linear steps, but make it coarser at larger radii so it looks better on a log plot
        plot_r = list(np.linspace(0.1, 1, 5))  #0.2kpc steps
        plot_r += list(np.linspace(1.5,10,17)) #0.5kpc steps
        plot_r += list(np.linspace(11,100,89)) #1.0kpc steps
        plot_r = np.array(plot_r)

        #potential, dpotential = get_potential(cat, plot_r, gal_pos)
        param = get_param(cat, plot_r, gal_pos, plot_param)

        linestyle = '-'
        if 'isolated' not in run:
            #param *= cat.time**1
            pass

        print(param[0])
        ax.plot(plot_r, param, linestyle=linestyle, label=names[i])


    ax.legend()
    ax.set_xlabel("Radius (kpc)")
    ax.set_ylabel(plot_param)
    ax.set_xscale("log")
    ax.set_yscale("log")

    fig.savefig(f"plots/{plot_param}.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()
