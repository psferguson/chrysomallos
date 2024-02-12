from lsst.daf.butler import Butler

__all__ = ["get_data_ids", "get_coadds", "get_coadd_dict"]


def get_data_ids(tract, patch, bands, skymap="hsc_rings_v1"):
    """
    Constructs a dictionary of data IDs for retrieving coadded images.

    Parameters:
    - tract: Tract number as an integer.
    - patch: Patch number within the tract.
    - bands: A list of strings representing the photometric bands.
    - skymap: Name of the skymap to use, defaulting to 'hsc_rings_v1'.

    Returns:
    - A dictionary keyed by band, each containing a data ID dictionary.
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
    Retrieves coadded images and related data from the Butler repository.

    Parameters:
    - butler: An instance of the lsst.daf.butler.Butler class.
    - data_id_dict: A dictionary of data IDs for fetching coadd data.
    - bands: A list of strings representing the photometric bands.

    Returns:
    - A dictionary with coadded image data, keyed by band.
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
    """
    Ensures a coadd dictionary is available, fetching data if necessary.

    Parameters:
    - coadd_dict: An existing coadd dictionary, or None to fetch new data.
    - config: Configuration dictionary specifying repository and data details.

    Returns:
    - A dictionary containing coadded image data for each specified band.
    """
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
