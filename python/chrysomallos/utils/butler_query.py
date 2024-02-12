from lsst.daf.butler import Butler

__all__ = ["get_data_ids", "get_coadds", "get_coadd_dict"]


def get_data_ids(tract, patch, bands, skymap="hsc_rings_v1"):
    """
    Get data IDs based on band, skymap, tract, and patch information.

    Parameters
    ----------
    tract : int, optional
        Tract number. Default is 9615.ß
    patch : int, optional
        Patch number. Default is 3.
    """
    data_id_dict = {}
    for band in bands:
        data_id_dict[band] = {
            "band": band,
            "skymap": skymap,
            "tract": tract,
            "patch": patch,
        }
    return data_id_dict


def get_coadds(butler, data_id_dict, bands):
    """
    Fetch coadd data based on the specified bands and data IDs.
    """
    # TODO: double check deepCoadd_calexp vs deepCoadd
    coadd_dict = {band: {} for band in bands}
    for band in bands:
        image = butler.get("deepCoadd_calexp", dataId=data_id_dict[band])
        coadd_dict[band]["image"] = image
        coadd_dict[band]["wcs"] = image.getWcs()
        coadd_dict[band]["bbox"] = image.getBBox()
        coadd_dict[band]["psf"] = image.getPsf()
        coadd_dict[band]["photo_calib"] = image.getPhotoCalib()
    return coadd_dict


def get_coadd_dict(coadd_dict, config):
    if coadd_dict:
        coadd_dict = coadd_dict
    else:
        butler = Butler(
            config["pipelines"]["repo"],
            collections=config["pipelines"]["input_collections"],
        )
        data_id_dict = get_data_ids(
            tract=config["pipelines"]["tract"],
            patch=config["pipelines"]["patch"],
            bands=config["pipelines"]["bands"],
            skymap=config["pipelines"]["skymap"],
        )
        # grab deepCoaddd
        coadd_dict = get_coadds(
            butler=butler, data_id_dict=data_id_dict, bands=config["pipelines"]["bands"]
        )
    return coadd_dict
