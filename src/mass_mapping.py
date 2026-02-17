import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from lenspack.image.inversion import ks93
from lenspack.utils import bin2d
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS

def bboxarea_and_realarea(image, pixel_scale):
    """"
    Assumes that image is:
      - rectangular 2D
      - contains nans for empty pixels
      - contains floats for other pixels
    """
    pixel_area = (pixel_scale / 60)**2 #arcmin^2
    area = np.sum(~np.isnan(image)) * pixel_area
    bbox = np.shape(image)[0] * np.shape(image)[1] * pixel_area
    return area, bbox

def project_ra_dec_to_cartesian(ra,dec):
    """
    Taken (ra, dec) and projects these coords to a cartesian
    (x, y) frame in arcsec relative to the median (ra,dec)
    """
    #define reference point
    dec0 = np.nanmedian(dec) * np.pi / 180.0
    ra0 = np.nanmedian(ra)

    x = ((ra0 - ra) * np.cos(dec0)) * 3600  # east minus west plus; cos correction
    y = (dec - np.mean(dec)) * 3600 # arcsec
    return x, y

def rotate_shears(g1, g2, angle):
    """"
    Rotate the input shear with specified angle
    angle in radians
    """
    g1_rotated = g1 * np.cos(2 * angle) - g2 * np.sin(2 * angle)
    g2_rotated = g1 * np.sin(2 * angle) + g2 * np.cos(2 * angle)
    return g1_rotated, g2_rotated

def get_convergence(g1, g2):
    """"
    input: g1, g2 : shear components in rectangular shape
    transform a shear catalog into convergence using
    Kaiser-Squires method
    output: kappa E- and B-modes
    """
    npad = 10 * g1.shape[0]
    # Zeropad
    g1map_pad = np.pad(g1, npad, constant_values=0)
    g2map_pad = np.pad(g2, npad, constant_values=0)
    kappaE_pad, kappaB_pad = ks93(g1map_pad, g2map_pad)
    # Remove zero pad
    kappaE = kappaE_pad[npad:-npad, npad:-npad]
    kappaB = kappaB_pad[npad:-npad, npad:-npad]
    return kappaE, kappaB

def bin_shear(drz_image_dir, shear_cat_dir, mean_gal_per_bin, params):
    """"
    input: HST drz image, shear catalog and required source galaxies per bin

    Generates shear components g1 and g2 binned in cartesian projected deltaRA and
    deltaDEC coordinates.

    outputs:
    - bounds: min, max coords source galaxies in shear cat
    - pos: coords of the bins in cartesian projected frame [deltaRA, deltaDEC]
    - ngal: number of source galaxies per bin
    - gamma_binned: binned shear components [g1, g2]
    """

    #Load drz image and shear
    image = fits.open(drz_image_dir)
    header = image[params['fits_extension']].header
    drz_image = image[params['fits_extension']].data
    shears = fits.getdata(shear_cat_dir)

    if len(shears) < 3:
        print(len(shears), "SOURCE GALAXIES IS NOT ENOUGH")
        return

    # ra, dec: coordinates of the source galaxy
    ra  = shears["ra"]
    dec = shears["dec"]
    # project (ra,dec) on tangent cartesian plane (x,y) in arcsec
    x, y = project_ra_dec_to_cartesian(ra,dec)

    # Shear components in pixel coordinates
    g1_image = shears["gamma1"]
    g2_image = shears["gamma2"]

    #Rotate shears to (ra, dec) coords
    rotation_angle = header[params['orientation_header']] * np.pi / 180 # rad
    g1, g2 = rotate_shears(g1_image, g2_image, rotation_angle)

    # Aspect ratio of the image
    x_min, x_max = np.percentile(x, [1, 99])
    y_min, y_max = np.percentile(y, [1, 99])
    aspect_ratio = (x_max - x_min) / (y_max - y_min)

    wcs = WCS(header)
    pixel_scale = np.sqrt((wcs.pixel_scale_matrix ** 2).sum(axis=0))[0] * 3600

    # Bounding box area and real area of the image
    bboxarea, area = bboxarea_and_realarea(drz_image, pixel_scale)
    bboxarea_to_area = bboxarea / area

    # Number of bins for final mass map
    bins = np.sqrt(len(x) / mean_gal_per_bin * bboxarea_to_area)
    nbin = (round(bins * aspect_ratio ** 0.5), round(bins / aspect_ratio ** 0.5))

    # Bin shear components based on galaxy position
    ra_map, dec_map = bin2d(x, y, v=(x, y), npix=nbin)
    g1map, g2map = bin2d(x, y, v=(g1, g2), npix=nbin)
    ngal = bin2d(x, y, v=None, npix=nbin)
    ngal = np.where(ngal == 0, np.nan, ngal)

    bounds = [x_min, x_max, y_min, y_max]
    pos = [ra_map.flatten(), dec_map]
    gamma_binned = [g1map, g2map]
    return bounds, pos, ngal, gamma_binned, area


def SNRmap(gamma_binned, smoothing=1):
    g1map, g2map = gamma_binned[0], gamma_binned[1]
    # Recover convergence via Kaiser-Squires inversion (Fourier Transform)
    kappaE, kappaB = get_convergence(g1map, g2map)

    # We apply a gaussian filter for better visual representation
    kappaE_smoothed = gaussian_filter(kappaE, smoothing)
    kappaB_smoothed = gaussian_filter(kappaB, smoothing)
    snr_map = kappaE_smoothed / np.std(kappaB_smoothed)
    return snr_map


def plot_ngal_gamma_snr(ngal,
                        pos,
                        gamma_binned,
                        snr_map,
                        bounds,
                        plot=True,
                        outfile='wlmap.png',
                        title='title'):


    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle(title)

    # Plot number of galaxies/bin
    im0 = axes[0].imshow(ngal, cmap="jet", origin='lower', extent=bounds)
    axes[0].set_title("Ngal per bin")
    axes[0].set_xlabel(r"$\Delta$RA [arcsec]")
    axes[0].set_ylabel(r"$\Delta$Dec [arcsec]")
    axes[0].set_aspect("equal")
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    # Plot shear field (vector field)
    ra_map, dec_map = pos[0], pos[1]
    X = ra_map.flatten()
    Y = dec_map.flatten()
    g1, g2 = gamma_binned[0], gamma_binned[1]
    U = g1.flatten()
    V = g2.flatten()

    angle = 0.5 * np.arctan2(V, U) * 180 / np.pi
    mask = ~np.isnan(U) & ~np.isnan(V)
    X, Y, U, V, angle = X[mask], Y[mask], U[mask], V[mask], angle[mask]
    axes[1].quiver(X, Y, U, V, angles=angle, pivot='middle', headaxislength=0,
                   headwidth=0, headlength=0, scale=0.005, scale_units='x', width=0.006)
    # axes[1].axis([np.min(y), np.max(y), np.min(x), np.max(x)])
    axes[1].set_title("Shear Field")
    axes[1].set_xlabel(r"$\Delta$RA [arcsec]")
    axes[1].set_ylabel(r"$\Delta$Dec [arcsec]")
    axes[1].set_aspect("equal")
    # Plot signal-to-noise map

    im2 = axes[2].imshow(snr_map,
                         extent=bounds,
                         cmap="jet", origin='lower', interpolation='quadric')

    axes[2].set_title("Signal-to-Noise Map")
    axes[2].set_xlabel(r"$\Delta$RA [arcsec]")
    axes[2].set_ylabel(r"$\Delta$Dec [arcsec]")
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    plt.tight_layout()

    #plt.savefig(outfile)
    #if plot:

    plt.show(block=True)

def generate_binned_lensing_map(drz_image_dir,
                                shear_cat_dir,
                                mean_gal_per_bin,
                                params):

    bounds, pos, ngal, gamma_binned, area = bin_shear(drz_image_dir,
                                                shear_cat_dir,
                                                mean_gal_per_bin,
                                                params)

    snr_map = SNRmap(gamma_binned, smoothing=1)

    total_ngals_density = np.nansum(ngal) / area

    plot_ngal_gamma_snr(ngal,
                        pos,
                        gamma_binned,
                        snr_map,
                        bounds,
                        plot=True,
                        outfile='wlmap.png',
                        title=r'Source galaxy number density: %.2f gal/arcmin^2' % total_ngals_density)
    return









