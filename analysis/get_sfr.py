import numpy as np
from readData.DataLoader import DataLoader
import units.springel_units as units  
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic

def main():

    path = "/home/j.rose/Projects/SMUGGLE/RUNs/CDM_du10/output-orange/"

    cat = DataLoader(path, snap_num=123, part_types=4, keys=['GFM_StellarFormationTime', 'Masses'])

    age = 13.7 - units.age_from_a(cat['GFM_StellarFormationTime'], H0=69.09, Om0=0.301712)
    cut = age > 0

    masses = cat['Masses']*1e10/cat.h

    print(np.min(age), np.max(age))

    np.save("plots/sfr", age[cut])
    np.save("plots/masses", masses[cut])

    sfr_bins = np.linspace(5.0, 7.5, 1000)
    mass_binned, bin_edges, bin_idx = binned_statistic(age, masses, 'sum', bins=sfr_bins)
    bin_half_width = (bin_edges[1]-bin_edges[0])/2
    bins = bin_edges[:-1] + bin_half_width
    sfr = mass_binned / (bin_half_width*2*1e9)

    fig, ax = plt.subplots()
    ax.plot(bins, sfr, label='Isolated')
    fig.savefig("plots/sfr.pdf")



    return

if __name__=="__main__":
    main()
