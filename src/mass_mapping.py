import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from lenspack.image.inversion import ks93
from lenspack.utils import bin2d
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS
from regions import Regions, PolygonSkyRegion, CircleSkyRegion, RectangleSkyRegion
from shapely.geometry import Polygon
from shapely.ops import unary_union
import astropy.units as u
import matplotlib as mpl
import re

def circle_to_polygon(circle, npoints=50):
    ra0 = circle.center.ra.deg
    dec0 = circle.center.dec.deg
    r = circle.radius.to(u.deg).value
    theta = np.linspace(0, 2*np.pi, npoints)
    x = ra0 + r * np.cos(theta) * np.cos(np.deg2rad(dec0))  # flat-sky
    y = dec0 + r * np.sin(theta)
    return Polygon(np.column_stack((x, y)))

def rectangle_to_polygon(rect):
    ra0 = rect.center.ra.deg
    dec0 = rect.center.dec.deg
    w = rect.width.to(u.deg).value
    h = rect.height.to(u.deg).value
    angle = rect.angle.to(u.rad).value  # rotation CCW from North

    # rectangle corners relative to center
    dx = np.array([-w/2, w/2, w/2, -w/2])
    dy = np.array([-h/2, -h/2, h/2, h/2])

    # rotate
    x = dx*np.cos(angle) - dy*np.sin(angle)
    y = dx*np.sin(angle) + dy*np.cos(angle)

    x += ra0
    y += dec0
    return Polygon(np.column_stack((x, y)))

def safe_read_ds9(mask_file):
    """
    Reads DS9 region file but filters invalid geometries before parsing.
    """

    valid_lines = []

    number_pattern = re.compile(r"-?\d+\.?\d*")

    for line in open(mask_file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        if line.startswith("circle"):
            nums = number_pattern.findall(line)
            if len(nums) >= 3:
                radius = float(nums[2])
                if radius <= 0:
                    print(f"Skipping invalid circle (radius<=0): {line}")
                    continue

        elif line.startswith("box"):
            nums = number_pattern.findall(line)
            if len(nums) >= 4:
                w = float(nums[2])
                h = float(nums[3])

                # fix or skip invalid boxes
                if w <= 0 or h <= 0:
                    print(f"Skipping invalid box (w/h<=0): {line}")
                    continue
        valid_lines.append(line)

    return Regions.parse("\n".join(valid_lines), format="ds9")

def get_masked_area(mask_file):
    """
    Reads a DS9 region file containing masks and calculates the total masked area in arcmin^2.
    """
    # Read regions from the DS9 region file
    #regions = Regions.read(mask_file, format="ds9")
    regions = safe_read_ds9(mask_file)
    polys = []
    for r in regions:
        if isinstance(r, PolygonSkyRegion):
            ra = r.vertices.ra.deg
            dec = r.vertices.dec.deg
            x = ra * np.cos(np.deg2rad(np.mean(dec)))
            y = dec
            polys.append(Polygon(np.column_stack((x, y))))
        elif isinstance(r, CircleSkyRegion):
            polys.append(circle_to_polygon(r))
        elif isinstance(r, RectangleSkyRegion):
            polys.append(rectangle_to_polygon(r))
        else:
            print(f"Warning: Skipping unsupported region type: {type(r)}")

    #Fix random errors in polygon creation
    fixed_polys = []
    for p in polys:
        if not p.is_valid:
            p = p.buffer(0)  # this often fixes self-intersections
        fixed_polys.append(p)

    # merge overlapping polygons
    merged = unary_union(fixed_polys)

    # area in deg^2, then convert to arcmin^2
    area_deg2 = merged.area
    area_arcmin2 = area_deg2 * 3600

    return area_arcmin2

def set_bw_dark_theme():
    """
    Pure black background + white/gray foreground (no color).
    Call once before creating figures.
    """
    plt.style.use("dark_background")

    mpl.rcParams.update(
        {
            # Surfaces
            "figure.facecolor": "#000000",
            "axes.facecolor": "#000000",
            "savefig.facecolor": "#000000",

            # Foreground/text
            "text.color": "#ffffff",
            "axes.labelcolor": "#ffffff",
            "axes.edgecolor": "#ffffff",
            "xtick.color": "#ffffff",
            "ytick.color": "#ffffff",

            # Lines/markers default (white)
            "lines.color": "#ffffff",
            "patch.edgecolor": "#ffffff",

            # Grid (subtle gray)
            "grid.color": "#666666",
            "grid.alpha": 0.35,
            "axes.grid": False,

            # Legend (B/W)
            "legend.frameon": True,
            "legend.facecolor": "#000000",
            "legend.edgecolor": "#ffffff",

            # Images/colormaps default to gray
            "image.cmap": "gray",
        }
    )
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

def project_ra_dec_to_cartesian(ra,dec, refcoord = None):
    """
    Taken (ra, dec) and projects these coords to a cartesian
    (x, y) frame in arcsec relative to the median (ra,dec)
    """
    #define reference point
    if refcoord is None:
        dec0 = np.nanmedian(dec)
        ra0 = np.nanmedian(ra)
    else:
        dec0 = refcoord[1]
        ra0 = refcoord[0]
    x = ((ra0 - ra) * np.cos(np.deg2rad(dec0))) * 3600  # east minus west plus; cos correction
    y = (dec - dec0) * 3600 # arcsec
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

def bin_shear(drz_image_dir, shear_cat_dir, mean_gal_per_bin, params, refcoord =None):
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
    x, y = project_ra_dec_to_cartesian(ra,dec, refcoord)

    # Shear components in pixel coordinates
    g1_image = shears["gamma1"]
    g2_image = shears["gamma2"]

    #Rotate shears to (ra, dec) coords
    if params['telescope'] == 'HST':
        rotation_angle = header[params['orientation_header']] * np.pi / 180 # rad
    elif params['telescope'] == 'JWST':
        rotation_angle = (header[params['orientation_header']]) * np.pi / 180
    g1, g2 = rotate_shears(g1_image, g2_image, rotation_angle)

    # Aspect ratio of the image
    x_min, x_max = np.percentile(x, [1, 99])
    y_min, y_max = np.percentile(y, [1, 99])
    aspect_ratio = (x_max - x_min) / (y_max - y_min)

    wcs = WCS(header)
    pixel_scale = np.sqrt((wcs.pixel_scale_matrix ** 2).sum(axis=0))[0] * 3600

    # Bounding box area and real area of the image
    area, bboxarea = bboxarea_and_realarea(drz_image, pixel_scale)
    bboxarea_to_area =  bboxarea / area

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


    #fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    set_bw_dark_theme()
    fig = plt.figure(figsize=(12, 4), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1])  # all equal
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    axes = [ax0, ax1, ax2]
    fig.suptitle(title)

    # Plot number of galaxies/bin
    im0 = axes[0].imshow(ngal, cmap="jet", origin='lower', extent=bounds)
    axes[0].set_title("Ngal per bin")
    axes[0].set_xlabel(r"$\Delta$RA [arcsec]")
    axes[0].set_ylabel(r"$\Delta$Dec [arcsec]")
    axes[0].set_aspect("equal")
    plt.colorbar(im0, ax=axes[0], fraction=0.049, pad=0.04)

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
    axes[1].quiver(X, Y, U, V, angles=angle, pivot='middle', headaxislength=0,color='white',
                   headwidth=0, headlength=0, scale=0.005, scale_units='x', width=0.006)
    # axes[1].axis([np.min(y), np.max(y), np.min(x), np.max(x)])
    axes[1].set_title("Shear Field")
    axes[1].set_xlabel(r"$\Delta$RA [arcsec]")
    axes[1].set_ylabel(r"$\Delta$Dec [arcsec]")
    axes[1].set_aspect("equal")
    #plt.colorbar(plt.cm.ScalarMappable(cmap="jet"), ax=axes[1], fraction=0.046, pad=0.04)
    # Plot signal-to-noise map

    im2 = axes[2].imshow(snr_map,
                         extent=bounds,
                         cmap="jet", origin='lower', interpolation='quadric')

    axes[2].set_title("Signal-to-Noise Map")
    axes[2].set_xlabel(r"$\Delta$RA [arcsec]")
    axes[2].set_ylabel(r"$\Delta$Dec [arcsec]")
    plt.colorbar(im2, ax=axes[2], fraction=0.049, pad=0.04)

    #Set lims
    for ax in axes:
        ax.set_xlim(bounds[0], bounds[1])
        ax.set_ylim(bounds[2], bounds[3])

    plt.show(block=True)

def eps2xi(eps_1, eps_2):
    """
    Convert eps to xi (eq 4.11 Bartellman & Schneider 2001)
    """
    eps_norm = np.sqrt(eps_1**2 + eps_2**2)
    xi_1 = 2*eps_1/(1 + eps_norm**2)
    xi_2 = 2*eps_2/(1 + eps_norm**2)
    return xi_1, xi_2

def generate_binned_lensing_map(drz_image_dir,
                                shear_cat_dir,
                                mean_gal_per_bin,
                                params):

    bounds, pos, ngal, gamma_binned, area = bin_shear(drz_image_dir,
                                                shear_cat_dir,
                                                mean_gal_per_bin,
                                                params)
    gamma_binned = [gamma_binned[0], gamma_binned[1]]
    snr_map = SNRmap(gamma_binned, smoothing=1)

    if params['mask_file'] is not None and os.path.isfile(params['mask_file']):
        masked_area = get_masked_area(params['mask_file'])
    else:
        masked_area = 0

    eff_area = area - masked_area
    print(f"Total area of the image: {area:.2f} arcmin^2 | Masked area: {masked_area:.2f} arcmin^2")
    total_ngals_density = np.nansum(ngal) / eff_area

    shears = fits.getdata(shear_cat_dir)
    eps1, eps2 = eps2xi(shears["e1"], shears["e2"])
    rms = 0.5 * (np.sqrt(np.mean(eps1 ** 2)) + np.sqrt(np.mean(eps2 ** 2)))

    plot_ngal_gamma_snr(ngal,
                        pos,
                        gamma_binned,
                        snr_map,
                        bounds,
                        plot=True,
                        outfile='wlmap.png',
                        title=f'Source galaxy number density: {total_ngals_density:.2f} gal/arcmin2 | RMS: {rms:.2f}')
    return









