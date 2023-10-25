import os
from starlink import DwarfConfig,CreateDwarfInjectionCatalog
import lsst.source.injection as si

default_config_dict={
    "repo" : "/repo/main",
    "collection" : 'HSC/runs/RC2/w_2023_32/DM-40356',
    "tract": 9615,
    "patch": 3,
    "bands" : ["g","r","i"],
    "dwarfs": []
    }
default_dwarf_dict={ 'id': 0,
                     'x_cen': 1000,
                     'y_cen': 1000,
                     'age': 10.0,
                     'feh': -2.0,
                     'mass': 500000.0,
                     'dist': 2.0,
                     'r_scale': 200,
                     'ellip': 0,
                     'theta': 0,
                     'n': 1,
                     'mag_limit': 27,
                     'mag_limit_band': 'LSST_g',
                     'random_seed': 69}


if __name__ == "__main__":
    x_cen_vals = [500, 1500, 2500, 3500]
    y_cen_vals = [500, 1500, 2500, 3500]
    mag_lim_vals = [27, 34, 28 , 29, 30, 31, 32, 33]

    i=0
    rnd=1
    outpath="./16_patch_configs/"
    if os.path.exists(outpath):
        print("outpath exists")
    else:
        os.mkdir(outpath)
        
    for mag_lim in mag_lim_vals:
        i += 1
        print(i)
        new_config_dict=default_config_dict.copy()
        new_config_dict["inject_cat_collection"] = f"u/pferguso/maglim_16_test/maglim_{mag_lim}_round_1"
        j=0
        dwarf_dicts=[]
        for x_cen in x_cen_vals:
            for y_cen in y_cen_vals:
                dwarf_dict=default_dwarf_dict.copy()
                dwarf_dict['id'] = j
                j += 1
                dwarf_dict['mag_limit']= mag_lim
                dwarf_dict['x_cen'] = x_cen
                dwarf_dict['y_cen'] = y_cen
                dwarf_dict['random_seed'] = 10 * j
                dwarf_dicts.append(dwarf_dict)
        new_config_dict["dwarfs"]=dwarf_dicts 
        ingest = False
        if ingest = True:
            creator=CreateDwarfInjectionCatalog(new_config_dict)
        
            creator.save_config(f"./16_patch_configs/231011_16_maglim_{mag_lim}_round_{rnd}")
            creator.run(ingest=True)
            
        repo="/repo/main"

        # butler = Butler(repo, collections=collection)
        # bands=["g"]
        # for band in bands:
        #     dataid = {'band': band, 'skymap': 'hsc_rings_v1', 'tract': 9615, 'patch': 3}
        #     coadd_dict = {}
        #     coadd_dict[band] = butler.get('deepCoadd_calexp', dataId=dataid)
        #     for mag_lim in mag_lim_vals:
        #         injection_catalog_collections = [f"u/pferguso/maglim_16_test/maglim_{mag_lim}_round_1"]
        #         injection_catalog_refs =  butler.registry.get("injection_catalog", 
        #                                                      collections=injection_catalog_collections,
        #                                                      where=f"band='{band}'")
        #         injection_catalogs = [butler.get(ref) for ref in injection_catalog_refs]
        #         inject_config = si.CoaddInjectConfig()
        #         inject_task = si.CoaddInjectTask(config=inject_config)
        #         inject_output = inject_task.run(
        #             injection_catalogs=injection_catalogs,
        #             input_exposure=coadd_dict[band].clone(),
        #             psf=coadd_dict[band].getPsf(),
        #             photo_calib=coadd_dict[band].getPhotoCalib(),
        #             wcs=coadd_dict[band].getWcs(),
        #         )
        # si_coadd_g_lim28 = inject_output_lim28.output_exposure
        # si_cat_out_g_lim28 = inject_output_lim28.output_catalog
