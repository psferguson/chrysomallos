import pandas as pd
import os
from glob import glob

import numpy as np
from tqdm import tqdm

from chrysomallos.injection import (
    CreateDwarfInjectionCatalog,
    DwarfParamSampler,
    PostageStampGenerator,
)
from chrysomallos.utils import Config
from chrysomallos.utils.annotations import get_anotation_box

# Limit the number of threads for NumPy, OpenBLAS, MKL, etc.
os.environ["OMP_NUM_THREADS"] = "1"  # OpenMP
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # OpenBLAS
os.environ["MKL_NUM_THREADS"] = "1"  # MKL
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # macOS Accelerate
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # NumExpr


def run_single_config(args):
    tract, patch, sampling_seed, config, patch_frame = args
    np.random.seed(sampling_seed)
    
    injection_bool = patch_frame['injection_bool'].values
    injection_frame = patch_frame[injection_bool] 
    single_config = config.copy()
    single_config['pipelines']['tract'] = tract
    single_config['pipelines']['patch'] = patch
    single_config['sampling']['n_dwarfs'] = len(injection_frame)
    config["sampling"]["generation_id"] = np.random.randint(0, 1000000)
    version = config["stamp"]["version"]
    config["stamp"]["title_format"] = f"madcash_stamp_v{version}"
    x_injects = injection_frame['x_cen'].values.copy()
    y_injects = injection_frame['y_cen'].values.copy()
    x_injects += np.random.randint(-injection_frame['cutout_size_x'].values//2, injection_frame['cutout_size_x'].values//2, len(injection_frame))
    y_injects += np.random.randint(-injection_frame['cutout_size_y'].values//2, injection_frame['cutout_size_y'].values//2, len(injection_frame))
    
    single_config['stamp']['stamp_x_cen'] = injection_frame['x_cen'].values
    single_config['stamp']['stamp_y_cen'] = injection_frame['y_cen'].values
    single_config['stamp']['stamp_indexes'] = injection_frame['stamp_index'].values
    single_config['sampling']['params']['x_cen'] = list(x_injects)
    single_config['sampling']['params']['y_cen'] = list(y_injects)
    single_config['sampling']['output_file'] = f'{tract}_{patch}_dwarf_injection_catalog.csv'
    if injection_bool.sum() > 0:   
        sampler = DwarfParamSampler(single_config)
        dwarf_params_frame, coadd_dict = sampler.run()

        creator = CreateDwarfInjectionCatalog(single_config, dwarf_params_frame, coadd_dict)
        print(f'{tract},{patch}: Creating dwarf injection catalog')
        catalogs, coadd_dict = creator.run(
            ingest=False,
            multiproc=False,
            crop_to_stamp=True,
        )

        postage_stamp_generator = PostageStampGenerator(
            config=single_config,
            dwarf_params_frame=dwarf_params_frame,
            dwarf_catalogs=catalogs,
            coadd_dict=coadd_dict,
        )
        print(f'{tract},{patch}: Generating postage stamps')
        postage_stamp_generator.run()
    else:
        postage_stamp_generator = PostageStampGenerator(
            config=single_config,
            dwarf_params_frame=None,
            dwarf_catalogs=None,
            coadd_dict=None,
        )
    empty_frame = patch_frame[~injection_bool]
    # check if empty_frame stamps are already generated
    # empty_frame_bool_dict={}
    # for band in single_config["pipelines"]["bands"]:
    #     empty_frame_bool_dict[band] = []
    #     for i in range(len(empty_frame)):
    #         stamp_dir = single_config["stamp"]["directory"] + "empty_stamps/" 
    #         os.makedirs(stamp_dir, exist_ok=True)
    #         stamp_title_prefix = single_config["stamp"]["title_format"]
    #         filename = stamp_dir + stamp_title_prefix
    #         filename += f"_{empty_frame['stamp_index'].values[i]}"
    #         filename += f"_band_{band}.fits"
    #         empty_frame_bool_dict[band].append(os.path.exists(filename))
    
    # empty_frame_bool = empty_frame_bool_dict[single_config["pipelines"]["bands"][0]] 
    # for band in single_config["pipelines"]["bands"][1:]:
    #     empty_frame_bool = np.logical_and(empty_frame_bool, empty_frame_bool_dict[band])
    # empty_frame= empty_frame[~np.array(empty_frame_bool)]
    # import pdb; pdb.set_trace()
    if len(empty_frame) > 0:
        postage_stamp_generator.generate_empty_stamps(
            n_stamps=len(empty_frame),
            stamp_x_cens=empty_frame['x_cen'].values,
            stamp_y_cens=empty_frame['y_cen'].values,
            stamp_indexes=empty_frame['stamp_index'].values,
        )
        
        print(f"{tract},{patch}: Generating {len(empty_frame)} empty stamps")

def run_configs(manifest_frame, config_dict, multiproc=False):
    
    gen_args = []
    config = Config(config_dict)

    tracts = manifest_frame['tract'].unique()

    for tract in tracts:
        patches = manifest_frame[manifest_frame['tract'] == tract]['patch'].unique()
        for patch in patches:
            patch_frame = manifest_frame[
                (manifest_frame['tract'] == tract) &\
                (manifest_frame['patch'] == patch)
            ]
            random_seed = tract * patch 
            gen_args.append(
                (tract, patch, random_seed, config, patch_frame)
            )
    if multiproc:
        from multiprocessing import Pool

        multiproc = int(multiproc)
        processes = multiproc if multiproc > 0 else None
        with Pool(processes, maxtasksperchild=1) as p:
            for _ in tqdm(
                p.imap_unordered(run_single_config, gen_args), 
                total=len(gen_args), 
                desc="Processing (multiprocess)"
            ):
                pass
        
    else:
        for args in tqdm(gen_args, desc="Processing (single process)"):
            run_single_config(args)


if __name__ == "__main__":
    main_dir = '/Volumes/gimli/hsc_pdr3/pferguson/dwarf_finder/'
    config_dict = "./madcash_v2_config.yaml"
    
    manifest_frame_path = main_dir + 'data/madcash/run_v2_stamp_centers.csv'
    manifest_frame = pd.read_csv(manifest_frame_path)
    run_configs(manifest_frame, config_dict, multiproc=15)
