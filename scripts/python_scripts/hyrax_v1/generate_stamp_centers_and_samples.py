import os
import pandas as pd
from math import ceil
from lsst.daf.butler import Butler
import numpy as np
import lsst.geom as geom
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

def get_stamp_centers(tract, patch, butler, cutout_size_x=600, cutout_size_y=600):
    pix_x_center = ceil(cutout_size_x / 2)
    pix_y_center = ceil(cutout_size_y / 2)
    
    calexp = butler.get("deepCoadd_calexp", dataId={'tract': tract, 
                                                    'patch': patch,
                                                    'band': 'r',
                                                    'skymap': 'hsc_rings_v1'}, 
                                                    collections=['pdr3_deepCoadd'])
    bbox = calexp.getBBox()
    wcs = calexp.getWcs()

    x_cens = []
    y_cens = []
    for x_cen in range(bbox.beginX + pix_x_center, bbox.endX, cutout_size_x):
        for y_cen in range(bbox.beginY + pix_y_center, bbox.endY, cutout_size_y):
            if bbox.endX - x_cen < pix_x_center:
                x_cen = bbox.endX - pix_x_center
            if bbox.endY - y_cen < pix_y_center:
                y_cen = bbox.endY - pix_x_center
            x_cens.append(x_cen - bbox.beginX)
            y_cens.append(y_cen - bbox.beginY)
    x_cens = np.asarray(x_cens)
    y_cens = np.asarray(y_cens)
    if (x_cens < 0).any() or (y_cens < 0).any():
        import pdb; pdb.set_trace()
        raise ValueError("x_cens and y_cens must be positive")
    # for each x_cen, y_cen, get the ra and dec
    ras = []
    decs = []
    for x_cen, y_cen in zip(x_cens, y_cens):
        pixel = geom.Point2D(x_cen, y_cen)
        sky = wcs.pixelToSky(pixel)
        ras.append(sky.getRa().asDegrees())
        decs.append(sky.getDec().asDegrees())
    
    injection_bool = np.random.uniform(0,1, len(x_cens)) < 0.5
    generation_seed = np.random.randint(0, 1000)
    inds = [f'{tract}_{patch}_{generation_seed:04d}_{i:03d}' for i in range(len(x_cens))]
    
    df = pd.DataFrame({
        'stamp_index': inds,
        'tract': [tract] * len(x_cens),
        'patch': [patch] * len(x_cens),
        'x_cen': x_cens,
        'y_cen': y_cens,
        'ra': ras,
        'dec': decs,
        'cutout_size_x': [cutout_size_x] * len(x_cens),
        'cutout_size_y': [cutout_size_y] * len(x_cens),
        'injection_bool': injection_bool,
    })
    return df

def process_patch(args):
    """Helper function to process a single patch."""
    tract, patch= args
    repo = "/Volumes/gimli/hsc_pdr3/repo/hsc_pdr3/"
    butler = Butler(repo)
    return get_stamp_centers(tract, patch, butler, cutout_size_x=600, cutout_size_y=600)

if __name__ == "__main__":
    
    repo = "/Volumes/gimli/hsc_pdr3/repo/hsc_pdr3/"
    butler = Butler(repo)
    registry = butler.registry

    coadd_refs = sorted(registry.queryDatasets('deepCoadd_calexp', collections = ['pdr3_deepCoadd']))
    skymap = butler.get("skyMap", dataId={'skymap':'hsc_rings_v1'}, collections=['skymaps'])
    # Extract tracts and patches
    tracts_patches = set([(i.dataId['tract'], i.dataId['patch']) for i in coadd_refs])
    tracts, patches = zip(*tracts_patches)

    ras, decs= [],[]
    for tract, patch in tracts_patches:
        tract_info= skymap[tract]
        wcs = tract_info.getWcs()
        patch_info = tract_info.getPatchInfo(patch)
        bbox = patch_info.getOuterBBox()
        center_pixel = geom.Point2D(bbox.getCenterX(), bbox.getCenterY())
        center_sky = wcs.pixelToSky(center_pixel)
        ras.append(center_sky.getRa().asDegrees())
        decs.append(center_sky.getDec().asDegrees())
    # Create a DataFrame
    patch_df = pd.DataFrame({
        'tract': tracts,
        'patch': patches,
        'ra': ras,
        'dec': decs
    })
    patch_df.to_csv('/Users/pferguson/projects/dwarf_finder/data/tract_patch_centers.csv', index=False)
    #create_injection_df
    injection_dfs = []
    args = [(tract, patch) for i, (tract, patch) in enumerate(tracts_patches)]
    with ProcessPoolExecutor(max_workers=15) as executor:
        results = list(tqdm(executor.map(process_patch, args), total=len(args), desc="Processing patches"))
        injection_dfs.extend(results)
    
    injection_df = pd.concat(injection_dfs, ignore_index=True)
    injection_df.to_csv('/Users/pferguson/projects/dwarf_finder/data/run_v1_stamp_centers.csv', index=False)

