import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
import shared_data as sd
from analysis.analyze import get_rotate_data
import h5py
import units.springel_units as units  

def save_pids(path, run, snap, name, part_types):
    
    keys = ['GroupPos', 'Group_R_Crit200', 'Coordinates', 'ParticleIDs']
    cat = DataLoader(path+run+'/output-orange/', snap_num=snap, part_types=part_types, keys=keys)

    gal_pos = cat['GroupPos'][0] 
    gal_r = cat['Group_R_Crit200'][0]
    print(gal_r * cat.time)
    for pt in part_types:
        coords = cat[f'PartType{pt}/Coordinates'] - gal_pos

        r = np.sqrt(np.sum(np.square(coords), axis=1))
        rcut = r < gal_r*2
        r2cut = r < gal_r*4

        np.save(f'plots/np_arrays/r200_{name}_{pt}_ids', cat[f'PartType{pt}/ParticleIDs'][rcut])
        if pt == 0:
            np.save(f'plots/np_arrays/2r200_{name}_{pt}_ids', cat[f'PartType{pt}/ParticleIDs'][r2cut])

    return

def make_ic_file(path, run, snap, name, part_types):

    with h5py.File(f"/home/j.rose//Projects/SMUGGLE/RUNs/CDM_du10/output-orange/snapdir_{snap:03}/snap_{snap:03}.0.hdf5", "r") as ofile:
        gas_keys = list(ofile['PartType0'].keys())
        dm_keys = list(ofile['PartType1'].keys())
        star_keys = list(ofile['PartType4'].keys())

    keys = list(set(gas_keys + dm_keys + star_keys))
    cat = DataLoader(path+run+'/output-orange/', snap_num=snap, part_types=[0,1,4], keys=keys)

    keys = {0:gas_keys, 1:dm_keys, 4:star_keys}

    with h5py.File(f"plots/ICs/{name}.hdf5", 'w') as ofile:

        for pt in part_types:
            group = ofile.create_group(f"PartType{pt}")

            ids = np.load(f"plots/np_arrays/r200_{name}_{pt}_ids.npy")
            
            cat_pt = pt
            if pt in [2, 3]:
                cat_pt = 4

            icut = np.in1d(cat[f'PartType{cat_pt}/ParticleIDs'], ids)
            if pt == 0:
                ids2 = np.load(f"plots/np_arrays/2r200_{name}_0_ids.npy")
                icut = np.in1d(cat['PartType0/ParticleIDs'], ids2)


            for key in keys[cat_pt]:
                data = cat[f'PartType{cat_pt}/{key}']
                if key == 'Coordinates':
                    data *= cat.time
                if key == 'Velocities':
                    data *= np.sqrt(cat.time)
                if key == 'Density':
                    data /= cat.time**3
                if key == 'Pressure':
                    data /= cat.time**3
                if key == 'GFM_StellarFormationTime':
                    data = 13.7 - units.age_from_a(data, H0=69.09, Om0=0.301712)
                if key == 'Potential':
                    data /= cat.time
                #if key == 'Softenings': 
                #    data *= cat.time #check this
    
                group.create_dataset(key, data=data[icut])


    return

def fix_header(path, run, snap, name, part_types):

    cat = DataLoader(path+run+'/output-orange/', snap_num=snap, part_types=-1)

    with h5py.File(f"plots/ICs/{name}.hdf5", 'r+') as ofile:

        group = ofile.create_group("Header")
        for key in cat.pheader:
            group.attrs[key] = cat.pheader[key]

        group.attrs['MassTable'] = cat.masstable
        group.attrs['NumFilesPerSnapshot'] = 1

        nparts = [len(ofile[f'PartType{pt}/Coordinates']) for pt in part_types]
        nparts_full = [0]*6
        for i,pt in enumerate(part_types):
            nparts_full[pt] = nparts[i]

        group.attrs['NumPart_ThisFile'] = np.array(nparts_full)
        group.attrs['NumPart_Total'] = np.array(nparts_full)

    return

def fix_cv(path, run, snap, name, part_types):

    keys = ['GroupPos', 'GroupVel'] 
    cat = DataLoader(path+run+'/output-orange/', snap_num=snap, part_types=-1, keys=keys)

    group_pos = cat['GroupPos'][0] * cat.time
    print(group_pos, cat.time)

    with h5py.File(f"plots/ICs/{name}.hdf5", 'r+') as ofile:

        size = 650
        for pt in part_types:
            coords = np.array(ofile[f'PartType{pt}/Coordinates']) - group_pos + size/2
            del ofile[f'PartType{pt}/Coordinates']  
            ofile[f'PartType{pt}'].create_dataset('Coordinates', data=coords)

        ofile['Header'].attrs['BoxSize'] = size
        print("Boxsize", size)

    return

def taper_outer_gas(path, run, snap, name, part_types):

    ids = np.load(f"plots/np_arrays/r200_0_ids.npy")
    ids2 = np.load(f"plots/np_arrays/2r200_0_ids.npy")

    with h5py.File(f"plots/ICs/{name}.hdf5", 'r+') as ofile:
        icut2 = np.in1d(np.array(ofile['PartType0/ParticleIDs']), ids2)
        icut1 = np.in1d(np.array(ofile['PartType0/ParticleIDs']), ids)
        icut = icut2 ^ icut1

        masses = np.array(ofile['PartType0/Masses'])
        temp = np.array(ofile['PartType0/InternalEnergy'])
        coords = np.array(ofile['PartType0/Coordinates'])[icut] - float(np.array(ofile['Header'].attrs['BoxSize']))/2
        
        r = np.sqrt(np.sum(np.square(coords), axis=1))
        rmi = 40
        rma = np.max(r)*1.3

        masses[icut] = masses[icut] * (rma - r)/(rma - rmi)
        temp[icut] = temp[icut] * (rma - r)/(rma - rmi) 

        masses[masses <= 0] = 1e-10
        temp[temp <= 0] = 1e-10

        del ofile['PartType0/InternalEnergy']
        del ofile['PartType0/Masses']

        ofile['PartType0'].create_dataset('InternalEnergy', data=temp)
        ofile['PartType0'].create_dataset('Masses', data=masses)

    return

def main():

    path='/home/j.rose/Projects/SMUGGLE/RUNs/'
    run = 'CDM_du10'
    snap = 75
    #name = 'du10'
    #name = 'du10_z2'
    #name = 'du10_z.4'
    name = 'du10_large2'

    save_pids(path, run, snap, name, [0,1,4])
    make_ic_file(path, run, snap, name, [0,1,4]) 
    fix_header(path, run, snap, name, [0,1,4]) 
    fix_cv(path, run, snap, name, [0,1,4]) 
    #taper_outer_gas(path, run, snap, name, 0)    

    return

if __name__=="__main__":
    main()
