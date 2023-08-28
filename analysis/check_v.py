from readData.DataLoader import DataLoader
import numpy as np
import matplotlib.pyplot as plt
import sys

from analysis import analyze
import shared_data

def main():

    path = "/home/j.rose/gr/"
    part_types = [4,0,1]
    keys = ["Coordinates","Velocities"] 

    G = 4.3e-6 #kpc/Msun * (km/s)^2

    path = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    run = 'CDM_du10_s75_phys_correct_time/output-blue/'
    all_runs = [path+run]

    path = '/home/j.rose/Projects/SMUGGLE/RUNs/'
    run = 'CDM_du10/output-orange'
    all_runs += [path + run]

    snaps = [0, 75]
    names = ['Isolated', 'Zoom']

    fig, ax = shared_data.set_plot_params(ncols=3, nrows=2)
    for i,run in enumerate(all_runs):
        vels = dict()
        coords = dict()
        if 'isolated' in run:
            cat = DataLoader(run, part_types=part_types, snap_num=snaps[i], keys=keys)
            for pt in part_types:
                vels[pt] = cat[f'PartType{pt}/Velocities'] 
                coords[pt] = cat[f'PartType{pt}/Coordinates'] 
        else:
            cat = DataLoader(run, part_types=part_types, snap_num=snaps[i], keys=keys+['GroupVel','GroupCM'], fof_idx=0)
            for pt in part_types:
                vels[pt] = cat[f'PartType{pt}/Velocities'] * np.sqrt(cat.time)
                coords[pt] = cat[f'PartType{pt}/Coordinates'] * cat.time

        if 'GroupCM' in cat:
            gal_pos = cat['GroupCM'] * cat.time 
            print(cat['GroupVel']/cat.time)
            print(cat['GroupCM']*cat.time)
        else:
            gal_pos = np.array([325]*3) 

        bins = np.linspace(-400, 400, 50)
        for j,pt in enumerate(part_types):
            box_cut = analyze.get_box_cut(coords[pt], gal_pos, 20)
            gal_vel = np.average(vels[pt][box_cut],axis=0)
            print(gal_vel)
            ax[i,j].hist(vels[pt][:,0] - gal_vel[0], bins=bins)
            ax[i,j].set_title(f"{names[i]} - PT{pt}")

    fig.savefig(f"plots/vhist.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()
