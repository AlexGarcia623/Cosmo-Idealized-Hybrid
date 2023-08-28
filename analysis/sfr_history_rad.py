from readData.DataLoader import DataLoader
import analysis.analyze as analyze
import shared_data
from scipy.interpolate import interp1d

import units.springel_units as units  

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.stats import binned_statistic

def main():

    part_types=[4]
    keys = ["GFM_StellarFormationTime", "GFM_InitialMass", "GroupSFR", 'GroupMassType', 'BirthPos']
    path = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    all_runs = ['CDM_du10_s75_phys_no_formation_time/output-blue/']
    #all_runs = ['CDM_du10_s61_phys_no_formation_time_no_fof/output-blue/']


    fig, ax = shared_data.set_plot_params(nrows=2, ncols=3)

    snaps = [171,127,127,127]
    #snaps = [snap_num]*len(names)

    names = ['CDM du10', 'CDM med', 'CDM med 2', 'CDM TNG']

    for i,run in enumerate(all_runs):

        cat = DataLoader(path+run, snaps[i], part_types=part_types, keys=keys, fof_idx=0) 

        if False:
            centers = []
            times = []
            for snap in range(61,snaps[i]+1):
                try:
                    cat = DataLoader(run, snap, keys='SubhaloPos',sub_idx=0)
                except:
                    try:
                        cat = DataLoader(run.split('-')[0] + '-orange', snap, keys='SubhaloPos',sub_idx=0)
                    except:
                        continue
                centers.append(cat['SubhaloPos'])
                times.append(cat.time)
            centers = np.array(centers)
            times = np.array(times)

            acut = (cat['GFM_StellarFormationTime'] > 0.4) & (cat['GFM_StellarFormationTime'] < np.max(times)-.01)

            #ages = units.age_from_a(cat['GFM_StellarFormationTime'][acut])  - units.age_from_a(1) # - units.age_from_a(cat.time)
            ages = cat.time
            masses = cat['GFM_InitialMass'][acut] * 1e10/cat.h

            #interpx_func = interp1d(times, centers[:,0])
            #interpy_func = interp1d(times, centers[:,1])
            #interpz_func = interp1d(times, centers[:,2])
            #adjusted_x = interpx_func(cat['GFM_StellarFormationTime'][acut])
            #adjusted_y = interpy_func(cat['GFM_StellarFormationTime'][acut])
            #adjusted_z = interpz_func(cat['GFM_StellarFormationTime'][acut])
            #adjusted_centers = np.array(list(zip(adjusted_x, adjusted_y, adjusted_z)))
            #pos = np.sqrt(np.sum(np.square(cat['BirthPos'][acut] - adjusted_centers), axis=1))
            #pos *= cat.time / cat.h

            cut = ~np.isnan(ages) #wind particles
            ages = ages[cut]
            masses = masses[cut]
            #pos = pos[cut]

            if False:
                cuts = [1,3,5,10,15,30]
                for j in range(len(cuts)): 
                    pcut = (pos < cuts[j])
                        
                    bins = np.linspace(0, 10, 150)

                    mass_binned, bin_edges, bin_idx = binned_statistic(ages[pcut], masses[pcut], 'sum', bins=bins)
                    
                    bin_half_width = (bin_edges[1]-bin_edges[0])/2
                    bins = bin_edges[:-1] + bin_half_width
                    sfr = mass_binned / (bin_half_width*2*1e9)
                    
                    ax[i//3,i%3].plot(bins, sfr)


        #ages = units.age_from_a(cat['GFM_StellarFormationTime'])  - units.age_from_a(1)
        ages = cat['GFM_StellarFormationTime']
        masses = cat['GFM_InitialMass'] * 1e10/cat.h
        cut = ages > 0 #wind particles
        ages = ages[cut]
        masses = masses[cut]
        bins = np.linspace(0, .2, 150)
        mass_binned, bin_edges, bin_idx = binned_statistic(ages, masses, 'sum', bins=bins)
        bin_half_width = (bin_edges[1]-bin_edges[0])/2
        bins = bin_edges[:-1] + bin_half_width
        sfr = mass_binned / (bin_half_width*2*1e9)
        p = ax[1,2].plot(bins, sfr, label=names[i])

        ax[i//3,i%3].set_title(f"{names[i]}")
        ax[i//3,i%3].set_xlabel("Look back time (Gyr)")
        ax[i//3,i%3].set_ylabel("SFR ($M_\odot$/yr)")
        #ax[i//3,i%3].set_xlim((14.5, 0))
        ax[i//3,i%3].set_ylim((0, 20))
    
    #for j in range(len(cuts)-1):
    #    ax[1,1].plot([0,0],[0,0], label=f"$<${cuts[j]}kpc")
    #ax[1,1].plot([0,0],[0,0], label=f"Total (in-situ)")
    #ax[1,1].legend()

    ax[1,2].set_xlabel("Look back time (Gyr)")
    ax[1,2].set_ylabel("SFR ($M_\odot$/yr)")
    ax[1,2].legend()
    #ax[1,2].set_xlim((14.5, 0))
    #ax[1,2].set_yscale("log")
    ax[1,2].set_ylim((0, 20))
    ax[1,2].set_title("Total (In \& Ex Situ)")

    fig.savefig("plots/sfr_history_rad.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()
