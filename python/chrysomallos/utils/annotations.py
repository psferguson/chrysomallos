import astropy.units as u
import lsst.geom as geom
import numpy as np

# def get_annotation_box():

# def create_annotation()


def get_anotation_box(wcs, bbox, x_cen, y_cen, r_scale, theta, ellip, scaling_factor=2):
    x0 = bbox.beginX
    y0 = bbox.beginY
    x_coord = x_cen + x0
    y_coord = y_cen + y0
    ra_cen, dec_cen = wcs.pixelToSky(x_coord, y_coord)

    # Calculate semi-major and semi-minor axes (units are delta arcsec)
    delta_a = scaling_factor * r_scale
    delta_b = delta_a * (1 - ellip)

    # Vertices of the ellipse before rotation
    vertices = np.array([[delta_a, 0], [-delta_a, 0], [0, delta_b], [0, -delta_b]])
    theta = np.radians(theta)
    # Rotation matrix
    rotation_matrix = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )

    # Rotate vertices
    rotated_vertices = np.dot(vertices, rotation_matrix.T)

    # Find bounding box
    min_ra, min_dec = np.min(rotated_vertices, axis=0)
    max_ra, max_dec = np.max(rotated_vertices, axis=0)

    # convert bounding box to pixel positions
    radec_min = geom.SpherePoint(
        ra_cen.asDegrees() + (min_ra * u.arcsec).to(u.deg).value,
        dec_cen.asDegrees() + (min_dec * u.arcsec).to(u.deg).value,
        geom.degrees,
    )
    radec_max = geom.SpherePoint(
        ra_cen.asDegrees() + (max_ra * u.arcsec).to(u.deg).value,
        dec_cen.asDegrees() + (max_dec * u.arcsec).to(u.deg).value,
        geom.degrees,
    )
    x_1, y_1 = wcs.skyToPixel(radec_min)
    x_2, y_2 = wcs.skyToPixel(radec_max)

    x_min = min(x_1, x_2)
    x_max = max(x_1, x_2)
    y_min = min(y_1, y_2)
    y_max = max(y_1, y_2)

    # need to confirm bounding box is within the image
    if x_min < bbox.beginX:
        x_min = bbox.beginX
    if x_max > bbox.endX:
        x_max = bbox.endX
    if y_min < bbox.beginY:
        y_min = bbox.beginY
    if y_max > bbox.endY:
        y_max = bbox.endY
    # but we want annotation in image coordinates so shift back
    x_min -= x0
    x_max -= x0
    y_min -= y0
    y_max -= y0

    return {
        "x_min": int(x_min),
        "x_max": int(x_max),
        "y_min": int(y_min),
        "y_max": int(y_max),
        "x_cen": x_cen,
        "y_cen": y_cen,
    }
