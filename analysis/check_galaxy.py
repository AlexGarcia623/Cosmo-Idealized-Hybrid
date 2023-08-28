import numpy as np
from readData.DataLoader import DataLoader 


def main():

    path = "/home/j.rose/Projects/SMUGGLE/isolated/RUNs/"
    snap = 0

    all_runs = ['MW_alex/output-blue']

    for i,run in enumerate(all_runs):

        cat = DataLoader(path + run, snap_num = snap, part_types=1, keys=['Coordinates'])

        print(np.min(cat['Coordinates'],axis=0))
        print(np.max(cat['Coordinates'],axis=0))
        print(np.average(cat['Coordinates'],axis=0))


    return

if __name__=="__main__":
    main()
