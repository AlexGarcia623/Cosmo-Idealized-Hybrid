import numpy as np
from readData.DataLoader import DataLoader as DL
import analysis.analyze as analyze
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def get_lambda(cat, npart, cuts, r200):

    G = 4.3e-6 #kpc/Msun * (km/s)^2

    coords = np.zeros((npart,3))
    vels = np.zeros((npart,3))
    masses = np.zeros(npart)
    start = 0
    for pt in [0,1,2,3]:
        num = np.sum(cuts[pt])
        coords[start:start+num] = cat[f'PartType{pt}/Coordinates'][cuts[pt]] - np.ones(3)*300
        vels[start:start+num]   = cat[f'PartType{pt}/Velocities'][cuts[pt]] - np.ones(3)*0
        masses[start:start+num] = cat[f'PartType{pt}/Masses'][cuts[pt]]*1e10
        start += num
    phi, theta = analyze.get_rotate_data(coords, vels, masses, face_on=True, get_pt=True, r_max=200)

    part_types = [0,1,2,3]
    mass = 0
    spec_lz = 0
    for pt in part_types:
        coords = cat[f'PartType{pt}/Coordinates'] - np.ones(3)*300
        vels = cat[f'PartType{pt}/Velocities'] - np.ones(3)*100
        masses = cat[f'PartType{pt}/Masses']*1e10
        coords, vels = analyze.get_rotate_data(coords, vels, masses, phi=phi, theta=theta)

        cut = np.sum(np.square(coords),axis=1) < r200**2
        mass += np.sum(masses[cut])
        spec_lz += np.sum((coords[cut,0]*vels[cut,1] - coords[cut,1]*vels[cut,0]) *masses[cut])/np.sum(masses[cut])
        #spec_lz += np.sum(coords[cut,0]*vels[cut,1] - coords[cut,1]*vels[cut,0]) 

    v200 = np.sqrt(G*mass/r200)
    lam = spec_lz / (np.sqrt(2) * v200 * r200)
    return lam

def get_r200(coords, mass, m200):

    r2 = np.sum(np.square(coords), axis=1)
    prev_m = 0
    for r in np.logspace(1,3,100):
        cut = r2 < r*r
        m = mass*np.sum(cut)
        if m > m200:
            break
        if m == prev_m:
            break
        prev_m = m
    return r

def lin(x, m, b):
    return m*x+b

def exp(x, a, b):
    return a*np.exp(b*x)

def calc_rs(coords, vels, masses):

    coords, vels = analyze.get_rotate_data(coords, vels, masses, face_on=True)
    part_r = np.sqrt(np.sum(np.square(coords), axis=1))

    dens_li = []
    all_r = np.logspace(-1,2,100)
    prev_r = 0
    for r in all_r:
        scut = (part_r < r) & (part_r > prev_r)
        mass = np.sum(masses[scut])

        out_area = np.pi*r**2
        in_area = np.pi*prev_r**2

        dens_li.append(mass / (out_area - in_area))
        prev_r = r

    r = np.array(all_r)
    d = np.array(dens_li) / 1e6

    fig, ax = plt.subplots()
    ax.plot(r,d)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim((1e-2,1e4))
    fig.savefig("/orange/paul.torrey/j.rose/rd.pdf")

    rcut = (r>5) & (np.exp(d) > 1e-1)
    guess = [-10, 8]
    (m,b), pcov = curve_fit(lin, r[rcut], d[rcut], guess)
    scale_len = -1/m

    return scale_len

def calc_height(coords, vels, masses):

    coords, vels = analyze.get_rotate_data(coords, vels, masses, edge_on=True)

    fig, ax = plt.subplots(figsize=(10,5))
    cut = analyze.get_box_cut(coords, [0,0,0], 30)
    xbins = np.linspace(-30,30,100)
    ybins = np.linspace(-5,5,100)
    ax.hist2d(coords[cut,1], coords[cut,0], norm=LogNorm(), bins=[xbins,ybins])
    fig.savefig("/orange/paul.torrey/j.rose/disc.pdf")


    return

def calc_conc(coords, masses):

    part_r = np.sqrt(np.sum(np.square(coords), axis=1))

    dens_li = []
    all_r = np.logspace(0,2,30)
    prev_r = 0
    for r in all_r:
        scut = (part_r < r) & (part_r > prev_r)
        mass = np.sum(masses[scut])

        out_area = np.pi*r**2
        in_area = np.pi*prev_r**2

        dens_li.append(mass / (out_area - in_area))
        prev_r = r   

    d = np.array(np.log10(dens_li))
    r = np.log10(all_r)
    di = np.diff(d)
    dd = di / np.diff(r)

    fig, ax = plt.subplots()
    ax.plot(r,d)
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.savefig("/orange/paul.torrey/j.rose/conc.pdf")

    idx = np.argmin(np.abs(dd+2))

    return np.power(10,r[idx])

def get_disc_cut(coords, vels, masses):


    coords, vels = analyze.get_rotate_data(coords, vels, masses, face_on=True)
    
    Lz = coords[:,0]*masses*vels[:,1] - coords[:,1]*masses*vels[:,0]
    Lmax = np.sqrt(np.sum(np.square(vels), axis=1)) * masses * np.sqrt(np.sum(np.square(coords), axis=1))

    Lmax[Lmax==0] = 0.1
    epsilon = Lz / Lmax

    return epsilon > .7


def main():

    path = "/home/j.rose/Projects/SMUGGLE/isolated/RUNs/"
    snap = 0
    #all_runs = ['MW_alex/output-blue/']
    all_runs = ['CDM_du10/output-blue/']

    #path = "/home/j.rose/Projects/SMUGGLE/RUNs/"
    #snap = 80
    #all_runs = ['CDM_du10/output-orange/']

    #disc height
    for run in all_runs:

        print(run, snap)

        if False:

            part_types = [0,1,4]

            cat = DL(path+run, part_types=part_types, snap_num=snap, keys=['Coordinates', 'Masses', 'Velocities', 'SubhaloMassType', 'SubhaloPos', 'SubhaloVel', 'Group_R_Crit200', 'Group_M_Crit200'], sub_idx=0, fof_idx=0)

            print(cat['SubhaloMassType'][5] * 1e10)

            coords = {pt:None for pt in part_types}
            vels = {pt:None for pt in part_types}
            masses = {pt:None for pt in part_types}
            for pt in part_types:
                coords[pt] = cat[f'PartType{pt}/Coordinates'] - cat['SubhaloPos']
                vels[pt] = cat[f'PartType{pt}/Velocities'] - cat['SubhaloVel']
                masses[pt] = cat[f'PartType{pt}/Masses'] * 1e10

            disc_cut = get_disc_cut(coords[4], vels[4], masses[4])
    
            coords[2] = coords[4][disc_cut]
            coords[3] = coords[4][~disc_cut]
            vels[2] = coords[4][disc_cut]
            vels[3] = coords[4][~disc_cut]
            masses[2] = masses[4][disc_cut]
            masses[3] = masses[4][~disc_cut]

            part_types = [0,1,2,3]
            nums = {pt:len(masses[pt]) for pt in part_types}

            h = cat.h
            redshift = cat.redshift

            #m200 = np.sum([np.sum(masses[pt]) for pt in part_types])
            m200 = cat['Group_M_Crit200'] * 1e10
            r200 = cat['Group_R_Crit200']

            npart = np.sum([np.sum(el) for el in nums])
            cuts = {pt:np.ones(len(masses[pt]),dtype=bool) for pt in part_types}
            #lam = get_lambda(cat, npart, cuts, r200)
            lam = 0

            #calculate radial scale length
            rs = calc_rs(coords[0], vels[0], masses[0])
           
            #calculate radial scale height
            height = calc_height(coords[0], vels[0], masses[0])

            #calculate concentration
            conc = calc_conc(coords[1], masses[1])

        else:

            cat = DL(path+run, part_types=[0,1,2,3], snap_num=snap, keys=['Coordinates', 'Masses', 'Velocities'])

            part_types = [0,1,2,3]

            center = np.average(cat[f'PartType1/Coordinates'], axis=0)
            bulk = np.average(cat[f'PartType1/Velocities'], axis=0)
            radius = 25

            #get particles inside the galaxy
            cuts = []
            for pt in part_types:
                pos = cat[f'PartType{pt}/Coordinates'] - center
                r2 = np.sum(np.square(pos), axis=1)
                if pt == 1:
                    cut = r2 < (radius * 5)**2
                else:
                    cut = r2 < radius**2
                cuts.append(cut)

            #get particle info
            masses = []
            pos = []
            vels = []
            nums = []
            for pt in part_types:
                masses.append(cat[f'PartType{pt}/Masses'][cuts[pt]]*1e10)
                pos.append(cat[f'PartType{pt}/Coordinates'][cuts[pt]] - center)
                vels.append(cat[f'PartType{pt}/Velocities'][cuts[pt]] - bulk)
                nums.append(np.sum(cuts[pt]))

            #read in certain values
            h = cat.h
            redshift = cat.redshift

            #calculate lambda
            m200 = 1.5e12
            m200 = np.sum([np.sum(masses[i]) for i in range(4)])
            r200 = get_r200(pos[1], masses[1][0], m200)
            print(r200)
            npart = np.sum(nums) 
            lam = get_lambda(cat, npart, cuts, r200)

            #calculate radial scale length
            rs = calc_rs(pos[0], vels[0], masses[0])
           
            #calculate radial scale height
            height = calc_height(pos[2], vels[2], masses[2])

            #calculate concentration
            conc = calc_conc(pos[1], masses[1])

        print(f"Masses: {np.sum(masses[0]):.2e}, {np.sum(masses[1]):.2e}, {np.sum(masses[2]):.2e}," \
                f"{np.sum(masses[3]):.2e}") #, {np.sum(masses[4]):.2e}")
        print(f"Particle Counts: {nums[0]:.2e}, {nums[1]:.2e}, {nums[2]:.2e}, {nums[3]:.2e}") #," \
        #        f"{nums[4]:.2e}")
        print("Hubble:", h)
        print("Redshift:", redshift)
        print(f"Virial Mass: {m200:.2e}")
        print(f"Lambda: {lam:.2f}")
        print(f"Radial Scale Length: {rs:.2f}")
        print(f"Concentration: {r200/conc:.1f}")

    return

if __name__=="__main__":
    main()
