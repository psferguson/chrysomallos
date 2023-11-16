import os
import sys
from math import ceil
from chrysomallos import DwarfConfig,CreateDwarfInjectionCatalog, sb_rh_to_mv, sb_mv_to_rh, mstar_from_absmag
import lsst.source.injection as si
import numpy as np
import matplotlib.pyplot as plt
from lsst.daf.butler import Butler, DimensionUniverse, DatasetType, CollectionType
import lsst.afw.display as afwDisplay
afwDisplay.setDefaultBackend('matplotlib')
from tqdm import tqdm
import astropy.units as u
import logging
logging.getLogger('chrysomallosLogger').setLevel(logging.CRITICAL)
logging.getLogger('lsst.coaddInjectTask').setLevel(logging.CRITICAL)

default_config_dict={
    "repo" : "/repo/main",
    "collection" : 'HSC/runs/RC2/w_2023_32/DM-40356',
    "tract": 9615,
    "patch": 3,
    "bands" : ["g"],
    "dwarfs": []
    }
default_dwarf_dict={ 'id': 0,
                     'x_cen': 1000,
                     'y_cen': 1000,
                     'age': 10.0,
                     'feh': -2.0,
                     'mass': 500000.0,
                     'dist': 4.0,
                     'r_scale': 200,
                     'ellip': 0,
                     'theta': 0,
                     'n': 1,
                     'mag_limit': 28,
                     'mag_limit_band': 'LSST_g'}

if __name__ == "__main__":
    repo="/repo/main"
    collection='HSC/runs/RC2/w_2023_32/DM-40356'
    butler = Butler(repo, collections=collection)
    registry = butler.registry

    mag_lim = 27
    x_cen=1500
    y_cen=2500
    dist = 2
    surface_brighness_vals = np.arange(23,27,1)
    mV_vals = np.arange(-4,-11,-2)

    new_config_dict=default_config_dict.copy()
    new_config_dict["inject_cat_collection"] = f"test"

    bands=["g"]
    for band in bands:
        dataid = {'band': band, 'skymap': 'hsc_rings_v1', 'tract': 9615, 'patch': 3}
        image_dict = {}
        image_dict[band] = butler.get('deepCoadd_calexp', dataId=dataid)

    wcs=image_dict["g"].getWcs()
    bbox=image_dict["g"].getBBox()

    catalogs={}
    j=0


    for sb in tqdm(surface_brighness_vals[:]):
        for mV in mV_vals:
            r_h = sb_mv_to_rh(sb,mV, distance=2e6)
            mass = mstar_from_absmag(mV)

            dwarf_dicts=[]
            new_config_dict=default_config_dict.copy()
            new_config_dict["inject_cat_collection"] = f"test"
            dwarf_dict=default_dwarf_dict.copy()
            dwarf_dict['id'] = j
            dwarf_dict["dist"]= dist
            j += 1
            dwarf_dict['x_cen'] = x_cen
            dwarf_dict['y_cen'] = y_cen
            dwarf_dicts.append(dwarf_dict)
            dwarf_dict["r_scale"]=r_h
            dwarf_dict["sb"]=sb
            dwarf_dict['mag_limit']=mag_lim

            dwarf_dict["m_v"]=mV


            dwarf_dict["mass"]=mass



            new_config_dict["dwarfs"]=dwarf_dicts
            creator=CreateDwarfInjectionCatalog(new_config_dict)
            tmp_catalog=creator.run(ingest=False,coadd_dict=image_dict)["g"]
            x_pix,y_pix = wcs.skyToPixelArray(tmp_catalog['ra'], tmp_catalog['dec'], degrees=True)
            sel = ((abs(x_pix-x_cen-bbox.beginX) < 400) &  (abs(y_pix-y_cen-bbox.beginY) < 400)) | (tmp_catalog["source_type"]=="Sersic")
            new_catalog_dict={"config":new_config_dict,"catalog":tmp_catalog[sel]}
            catalogs[j]=new_catalog_dict
    # plot samples
    x=[catalogs[i]["config"]["dwarfs"][0]["r_scale"] for i in catalogs.keys()]
    y=[catalogs[i]["config"]["dwarfs"][0]["m_v"] for i in catalogs.keys()]
    plt.figure()
    plt.scatter(x,y)
    plt.xscale('log')
    plt.ylim(5,-20)
    plt.xlabel("log r_scale")
    plt.ylabel("MV")
    plt.savefig("./parameters_space_MV_logr.png")



    inject_config = si.CoaddInjectConfig()
    inject_task = si.CoaddInjectTask(config=inject_config)

    # injection
    for i in tqdm(catalogs.keys()):
        inject_output = inject_task.run(
            injection_catalogs=[catalogs[i]["catalog"]],
            input_exposure=image_dict['g'].clone(),
            psf=image_dict['g'].getPsf(),
            photo_calib=image_dict['g'].getPhotoCalib(),
            wcs=image_dict['g'].getWcs(),
        )
        catalogs[i]["image"]=inject_output.output_exposure

    # create plot
    minx=bbox.beginX
    miny=bbox.beginY
    # plot_si_calexp.image.array = gaussian_filter(si_calexp.image.array, sigma=3)
    nrows=ceil((len(catalogs) ) / 4)
    fig, axs_arr = plt.subplots(nrows+1, 4, figsize=(12, 3 * nrows), dpi=150)
    axs=axs_arr.ravel()
    plt.sca(ax=axs[2])
    display0 = afwDisplay.Display(frame=fig)
    # display0.scale('asinh', min=-5/Q, max=25/Q, Q=Q)
    # display0.scale('linear', 'zscale')
    display0.scale('linear', min=-0.10, max=0.25)
    #display0.mtv(coadd_g.image)
    display0.mtv(image_dict['g'].image[minx+x_cen-300:minx+x_cen+300, y_cen-300:y_cen+300])
    plt.title('coadd image')

    for i in [0,1,3]:
        if i ==1:
            axs[i].set_title(f"injection dist = {dist} mpc")
        axs[i].set_axis_off()

    for i in tqdm(catalogs.keys()):
    #for i in tqdm(np.array([1,2,3,4])*4):
        scale=ceil(i / 4)
        plt.sca(axs[3+i])

        display0 = afwDisplay.Display(frame=fig)
        display0.scale('linear', min=-0.10, max=1/scale)

        display0.mtv(catalogs[i]["image"].image[minx+x_cen-300:minx+x_cen+300, y_cen-300:y_cen+300])
        title_str=f'log_rh={catalogs[i]["config"]["dwarfs"][0]["r_scale"]:0.0f} pc, '
        title_str+=f'Mv={catalogs[i]["config"]["dwarfs"][0]["m_v"]:0.1f}, '
        title_str+=f'\nsb={catalogs[i]["config"]["dwarfs"][0]["sb"]:0.1f},'
        plt.title(title_str)


    plt.tight_layout()
    plt.savefig(f"./sb_rh_test_{dist}_mpc.png")