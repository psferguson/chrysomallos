from tqdm import tqdm

class CatalogCreator:
    def __init__(self, butler, bands=["g"], default_config_dict=None, default_dwarf_dict=None):
        self.butler = butler
        self.bands = bands
        self.default_config_dict = default_config_dict or {}
        self.default_dwarf_dict = default_dwarf_dict or {}
        self.catalogs = {}
        self.image_dict = self._initialize_image_dict()

    def _initialize_image_dict(self):
        image_dict = {}
        for band in self.bands:
            dataid = {'band': band, 'skymap': 'hsc_rings_v1', 'tract': 9615, 'patch': 3}
            image_dict[band] = self.butler.get('deepCoadd_calexp', dataId=dataid)
        return image_dict

    def create_catalogs(self, surface_brighness_vals, mV_vals, dist, x_cen, y_cen, mag_lim):
        wcs = self.image_dict["g"].getWcs()
        bbox = self.image_dict["g"].getBBox()
        j = 0

        for sb in tqdm(surface_brighness_vals[:]):
            for mV in mV_vals:
                r_h = sb_mv_to_rh(sb, mV, distance=2e6)  # Assuming this function is defined elsewhere
                mass = mstar_from_absmag(mV)  # Assuming this function is defined elsewhere

                dwarf_dicts = []
                new_config_dict = self.default_config_dict.copy()
                new_config_dict["inject_cat_collection"] = "test"
                
                dwarf_dict = self.default_dwarf_dict.copy()
                dwarf_dict.update({
                    'id': j,
                    'dist': dist,
                    'x_cen': x_cen,
                    'y_cen': y_cen,
                    'r_scale': r_h,
                    'sb': sb,
                    'mag_limit': mag_lim,
                    'm_v': mV,
                    'mass': mass
                })

                dwarf_dicts.append(dwarf_dict)
                new_config_dict["dwarfs"] = dwarf_dicts 
                
                creator = CreateDwarfInjectionCatalog(new_config_dict)
                tmp_catalog = creator.run(ingest=False, coadd_dict=self.image_dict)["g"]
                
                x_pix, y_pix = wcs.skyToPixelArray(tmp_catalog['ra'], tmp_catalog['dec'], degrees=True)
                sel = ((abs(x_pix-x_cen-bbox.beginX) < 400) &  (abs(y_pix-y_cen-bbox.beginY) < 400)) | (tmp_catalog["source_type"]=="Sersic")
                new_catalog_dict = {"config": new_config_dict, "catalog": tmp_catalog[sel]}
                
                self.catalogs[j] = new_catalog_dict
                j += 1

    def get_catalogs(self):
        return self.catalogs
    
    
    
    
    
    def _create_catalog_for_params(self, params):
        sb, mV, wcs, bbox, dist, x_cen, y_cen, mag_lim = params

        r_h = sb_mv_to_rh(sb, mV, distance=2e6)
        mass = mstar_from_absmag(mV)

        dwarf_dicts = []
        new_config_dict = self.default_config_dict.copy()
        new_config_dict["inject_cat_collection"] = "test"
        
        dwarf_dict = self.default_dwarf_dict.copy()
        dwarf_dict.update({
            'id': self.j,
            'dist': dist,
            'x_cen': x_cen,
            'y_cen': y_cen,
            'r_scale': r_h,
            'sb': sb,
            'mag_limit': mag_lim,
            'm_v': mV,
            'mass': mass
        })

        dwarf_dicts.append(dwarf_dict)
        new_config_dict["dwarfs"] = dwarf_dicts 
        
        creator = CreateDwarfInjectionCatalog(new_config_dict)
        tmp_catalog = creator.run(ingest=False, coadd_dict=self.image_dict)["g"]
        
        x_pix, y_pix = wcs.skyToPixelArray(tmp_catalog['ra'], tmp_catalog['dec'], degrees=True)
        sel = ((abs(x_pix-x_cen-bbox.beginX) < 400) &  (abs(y_pix-y_cen-bbox.beginY) < 400)) | (tmp_catalog["source_type"]=="Sersic")
        new_catalog_dict = {"config": new_config_dict, "catalog": tmp_catalog[sel]}

        return (self.j, new_catalog_dict)

    def create_catalogs(self, surface_brighness_vals, mV_vals, dist, x_cen, y_cen, mag_lim):
        wcs = self.image_dict["g"].getWcs()
        bbox = self.image_dict["g"].getBBox()

        params = [(sb, mV, wcs, bbox, dist, x_cen, y_cen, mag_lim) 
                  for sb in surface_brighness_vals for mV in mV_vals]

        with multiprocessing.Pool() as pool:
            results = list(tqdm(pool.imap(self._create_catalog_for_params, params), total=len(params)))

        for j, catalog_dict in results:
            self.catalogs[j] = catalog_dict
            self.j += 1