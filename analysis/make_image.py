import matplotlib
matplotlib.use("agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from copy import copy

from readData.DataLoader import DataLoader 
from analysis import analyze

import visualization.contour_makepic as cmakepic
import util.calc_hsml as calc_hsml
import shared_data

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

    coords = copy(cat['Coordinates']) #*cat.time
    vels = copy(cat['Velocities'])
    metallicity = copy(cat['GFM_Metallicity'])
    sfr = copy(cat['StarFormationRate'])
    ea = copy(cat['ElectronAbundance'])
    ia = copy(cat['InternalEnergy'])
    
    masstab = cat.masstable[part_types[0]]       
    if masstab == 0:
        masses = cat['Masses']
    else:
        masses = np.array([masstab]*len(vels))    

    print(np.min(vels, axis=0), np.max(vels, axis=0))

    coords -= gal_pos
    vels -= gal_vel

    if 'isolated' not in cat.path:
        coords *= cat.time
        #vels *= np.sqrt(cat.time) 

    box_cut = analyze.get_box_cut(coords, np.zeros(3), maxfov)
    box_cut = np.ones(len(coords), dtype=bool)

    coords = coords[box_cut]
    vels = vels[box_cut]
    masses = masses[box_cut]
    metallicity = metallicity[box_cut]
    sfr = sfr[box_cut]
    ea = ea[box_cut]
    ia = ia[box_cut]

    print(np.min(vels,axis=0), np.max(vels,axis=0))

    #coords, vels = analyze.get_rotate_data(coords, vels, masses, face_on=face_on, edge_on=edge_on, r_min=1, r_max=7) 
    #phi, theta = analyze.get_rotate_data(coords, vels, masses, face_on=face_on, edge_on=edge_on, r_min=1, r_max=7, get_pt=True) 
    #print(phi, theta)
    #return
    coords, vels = analyze.get_rotate_data(coords, vels, masses, phi=1.07, theta=2.15) 

    print(np.min(vels,axis=0), np.max(vels,axis=0))

    if 'SubfindHsml' in cat.data:
        hsml = cat['SubfindHsml']
    else:
        hsml = calc_hsml.get_particle_hsml(coords[:,0], coords[:,1], coords[:,2])


    massmapx,image = cmakepic.simple_makepic(coords[:,0], coords[:,1],
            weights = masses, hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
    massmapy,image = cmakepic.simple_makepic(coords[:,0], coords[:,1],
            weights = metallicity, hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
    massmapz,image = cmakepic.simple_makepic(coords[:,0], coords[:,1],
            weights = np.abs(vels[:,0]), hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
    massmaps,image = cmakepic.simple_makepic(coords[:,0], coords[:,1],
            weights = sfr, hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
    massmape,image = cmakepic.simple_makepic(coords[:,0], coords[:,1],
            weights = ea, hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)
    massmapi,image = cmakepic.simple_makepic(coords[:,0], coords[:,1],
            weights = ia, hsml=hsml,xrange=xrange,yrange=yrange, pixels=pixels)

    return massmapx, massmapy, massmapz, massmaps, massmape, massmapi

def get_subplots(all_runs, do_faces):
    if do_faces:
        fig,ax = shared_data.set_plot_params(nrows=len(all_runs), ncols=6)
    else:
        fig,ax = shared_data.set_plot_params(nrows=len(all_runs), ncols=5)
    return fig, ax

def main():

    pixels = 512 #can be made larger #does this actually do anything?
    fov = 20 #kpc - diameter (not radius)
    do_faces = True

    ipath = '/home/j.rose/Projects/SMUGGLE/isolated/RUNs/'
    zpath = '/home/j.rose/Projects/SMUGGLE/RUNs/'
    keys = ["Coordinates", "Velocities", "Masses", "Potential", "GFM_Metallicity", "StarFormationRate", "ElectronAbundance", 'InternalEnergy'] 

    #all_runs = ['CDM_du10_s75_phys_no_formation_time/output-blue/']
    iso_runs = ['new_vel/output-blue/']
    all_runs = [ipath + run for run in iso_runs]

    zoom_runs = ['CDM_du10/output-orange/']
    #all_runs += [zpath + run for run in zoom_runs]

    snaps = [0, 75]
    names = ['Isolated', 'Zoom']

    fig, ax = get_subplots(all_runs, do_faces)
    if len(ax.shape) == 1:
        ax = np.array([ax])

    for i,run in enumerate(all_runs):
        print(run)
        ax_idx = 0
        snap_num = snaps[i]

        if 'isolated' not in run:
            cat = DataLoader(run, part_types=0, snap_num=snap_num, keys=keys+['GroupPos','GroupCM'], fof_idx=0) 
            gal_pos = cat['GroupPos']
            print(gal_pos)

        else:
            cat = DataLoader(run, part_types=[0], snap_num=snap_num, keys=keys)  
            gal_pos = np.average(cat['Coordinates'],axis=0)
            gal_vel = np.average(cat['Velocities'],axis=0)

        rads = [150,100,75,50,25]
        for rad in rads:
            box_cut = analyze.get_box_cut(cat['Coordinates'], gal_pos, rad)
            gal_pos = np.average(cat['Coordinates'][box_cut],axis=0,weights=cat['Masses'][box_cut])
        gal_vel = np.average(cat['Velocities'][box_cut],axis=0,weights=cat['Masses'][box_cut])

        print(gal_pos)

        if False:

            if do_faces:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[1], filt='', gal_pos=gal_pos, gal_vel=gal_vel)
            else:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=False, part_types=[1], filt='', gal_pos=gal_pos, gal_vel=gal_vel)

            im = ax[i,0].imshow(massmap, norm=LogNorm()) #vmin=1e6, vmax=1e9))
            ax[i,0].set_xticks([pixels/5*j for j in range(6)])
            ax[i,0].set_yticks([pixels/5*j for j in range(6)])
            ax[i,0].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,0].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])

            ax[i,0].set_title("DM")
            ax_idx += 1

        #cat = DataLoader(run, part_types=[0], snap_num=snap_num, keys=keys) 

        if True:

            if do_faces:
                massmaps = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[0], filt='', gal_pos=gal_pos, gal_vel=gal_vel)
            else:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=False, part_types=[0], filt='', gal_pos=gal_pos, gal_vel=gal_vel)


            vmin = []
            vmax = []
            titles = []
            for j,massmap in enumerate(massmaps):
                im = ax[i,ax_idx].imshow(massmap, norm=LogNorm())
                fig.colorbar(im, ax=ax[i,ax_idx])
                ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
                ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
                ax[i,ax_idx].set_xlabel('X [kpc]')
                ax[i,ax_idx].set_ylabel('Y [kpc]')

                ax[i,ax_idx].set_title(f"{names[i]} Gas Density")
                ax_idx += 1

        if False:
            cat = DataLoader(path+run, part_types=[0], snap_num=snap_num, keys=keys) 

            if do_faces:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=True, part_types=[0], filt='', gal_pos=gal_pos, gal_vel=gal_vel)

                ax[i,ax_idx].imshow(massmap, norm=LogNorm())
                ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
                ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])

                ax[i,ax_idx].set_title("Gas (Edge)")
                ax_idx += 1

        if False:

            cat = DataLoader(path+run, part_types=[2], snap_num=snap_num, keys=keys) 

            if do_faces:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[2], is_light=False, gal_pos=gal_pos, gal_vel=gal_vel)
            else:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=False, part_types=[2], is_light=False, gal_pos=gal_pos, gal_vel=gal_vel)

            ax[i,ax_idx].imshow(massmap, norm=LogNorm())
            ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_title("Disk Mass (Face)")   
            ax_idx += 1

            if do_faces:
                #cat = DataLoader(run, part_types=[4], snap_num=snap_num, keys=keys, sub_idx=-1, fof_idx=0) 
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=True, part_types=[3], is_light=False, gal_pos=gal_pos, gal_vel=gal_vel)

                ax[i,ax_idx].imshow(massmap, norm=LogNorm())
                ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
                ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
                ax[i,ax_idx].set_title("Disk Mass (Edge)")   
                ax_idx += 1


        if False:

            cat = DataLoader(path+run, part_types=[3], snap_num=snap_num, keys=keys) 

            if do_faces:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[3], is_light=False, gal_pos=gal_pos, gal_vel=gal_vel)
            else:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=False, part_types=[3], is_light=False, gal_pos=gal_pos, gal_vel=gal_vel)

            ax[i,ax_idx].imshow(massmap, norm=LogNorm())
            ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_title("Bulge Mass")   
            ax_idx += 1
        
        if False:

            cat = DataLoader(path+run, part_types=[4], snap_num=snap_num, keys=keys) 

            if do_faces:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=True, edge_on=False, part_types=[4], is_light=False)
            else:
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=False, part_types=[4], is_light=False)

            im = ax[i,ax_idx].imshow(massmap, norm=LogNorm(vmin=1e-6))
            ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
            ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
            ax[i,ax_idx].set_title("Stellar Mass (Face)")   
            #fig.colorbar(im, ax = ax[i,ax_idx])
            ax_idx += 1

            if do_faces:
                #cat = DataLoader(run, part_types=[4], snap_num=snap_num, keys=keys, sub_idx=-1, fof_idx=0) 
                massmap = get_massmap(cat, pixels=pixels, fov=fov, face_on=False, edge_on=True, part_types=[3], is_light=False)

                im = ax[i,ax_idx].imshow(massmap, norm=LogNorm(vmin=1e-6))
                ax[i,ax_idx].set_xticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_yticks([pixels/5*j for j in range(6)])
                ax[i,ax_idx].set_xticklabels([f'{fov/5*j:.0f}' for j in range(6)])
                ax[i,ax_idx].set_yticklabels([f'{fov/5*j:.0f}' for j in range(6)])
                ax[i,ax_idx].set_title("Stellar Mass (Edge)")   
                #fig.colorbar(im, ax = ax[i,ax_idx])

    fig.savefig(f"plots/light_map_{snap_num}.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()

