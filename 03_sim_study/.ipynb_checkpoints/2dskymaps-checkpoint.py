# imports
import argparse
import healpy as hp
import matplotlib as plt
from numpy import log, sqrt, cos, pi
import numpy as np
import math
import pylab

DEGREE = pi / 180.


def plot_map(skymap, title, label='', proj='C', dMin=None, dMax=None,
             filename=None, watermark='', half_map=False):
    """
    Plots a 2D skymap with a title. Can specify min, max of scale, filename
    and half map or full map.

    """
    # Colorway
    colormap = pylab.get_cmap("rainbow")
    notext = False

    if proj == 'C':
        rotation = (-180, 0, 0)
    else:
        rotation = (0, 0, 0)

    hp.mollview(skymap,
                fig=1,
                title=title,
                rot=rotation, coord=[proj],
                unit=label,
                notext=notext, cmap=colormap, min=dMin, max=dMax)

    fig = pylab.figure(1)
    for ax in fig.get_axes():
        if proj == 'C0':
            ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="x-large")
            ax.annotate("360$^\circ$", xy=(-1.95, 0.625), size="x-large")
        elif proj == 'C':
            ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="x-large")
            ax.annotate("360$^\circ$", xy=(-1.95, 0.625), size="x-large")
            if half_map:
                ax.annotate("0$^\circ$", xy=(1.8, -0.625), size="x-large")
                ax.annotate("360$^\circ$", xy=(-1.99, -0.625), size="x-large")

    hp.graticule()
    hp.projtext(95 * DEGREE + (rotation[0] - 180) * DEGREE,
                280 * DEGREE - rotation[1] * DEGREE,
                watermark,
                # coord=coords,
                color="grey",
                alpha=0.5,
                rotation=0,
                fontdict={"family": "sans-serif", "weight": "bold", "size": 42})

    if half_map:
        plt.ylim([-1, 0])

    if filename:
        fig.savefig(filename, transparent=True, dpi=100)

    plt.show()


def mask_map(map, dec_min, dec_max):
    """
    Masks map to cut data outside of specified min and max declination (inclusive).
    :param map: healpy map to be masked.
    :param dec_min: Minimum declination to include (inclusive)
    :param dec_max: Maximum declination to include (inclusive)
    :return: Masked map.
    """
    npix = len(map)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))

    thetaMin = (90 - dec_min) * DEGREE
    thetaMax = (90 - dec_max) * DEGREE
    thetaCut = (theta <= thetaMin) * (theta >= thetaMax)

    new_map = np.copy(map)  # not copy.deepcopy(map)?
    new_map[np.logical_not(thetaCut)] = hp.UNSEEN

    return new_map


def top_hat_smooth(orig_map, radius, average=False):
    """
    Applies Top Hat Smoothing algorithm to map with specified radius in degrees.

    :param orig_map: healpy map. Original unsmoothed map.
    :param radius: float, Radius of smoothing in degrees.
    :param average: Boolean. True averages pixels. False sums surrounding pixels.
    :return: healpy map smoothed
    """
    nside = hp.npix2nside(len(orig_map))
    pixel_list = []
    dsum = 0.0
    new_map = np.zeros(len(orig_map))

    # Loop over all pixels
    for ipix in range(0, len(orig_map)):
        # Added if statement - 01_30_2024
        if orig_map[ipix] == hp.UNSEEN:
            new_map[ipix] = hp.UNSEEN
        else:
            theta, phi = hp.pix2ang(nside, ipix)
            center = hp.dir2vec(theta, phi)

            # Grab pixels with given angular radius of ith pixel
            pixel_list = hp.query_disc(nside, center, radius)

            # Sum up pixels in the angular cap; set value of ith pixel to this sum
            dsum = 0.
            nsum = 0.
            for jpix in pixel_list:
                if orig_map[jpix] != hp.UNSEEN:
                    dsum += orig_map[jpix]
                    nsum += 1.0

            if average and nsum > 0:
                new_map[ipix] = dsum / nsum  # set to average
            else:
                new_map[ipix] = dsum

    return new_map


def relative_intensity(map, title, dec_min, dec_max, smoothing_angle, pngname,
                       min=None, max=None):
    """
    Plots 2D half map of relative intensity, masking and then smoothing, and
    saves it as a png file in a specified output directory.
    Scale can optionally be specified.

    :param map: String path to fits file with  relative intensity.
    :param title: String title of plot.
    :param dec_min: float minimum declination for masking.
    :param dec_max: float maximum declination for masking.
    :param smoothing_angle: float radius for smoothing.
    :param pngname: String path of png file that will be saved
    :param min: float minimum of scale. Optional, default is none.
    :param max: float maximum of scale. Optional, default is none.

    :return: Saves png file of skymap.
    """
    masked_map = mask_map(map, dec_min, dec_max)
    smoothed_map = top_hat_smooth(masked_map, smoothing_angle * DEGREE,
                                  average=True)

    # Plot with specified range and save
    plot_map(smoothed_map, title, label='Rel. Int. ', half_map=True, dMin=min,
             dMax=max, filename=pngname)


def li_ma_significance(rel_int, variance, title, dec_min, dec_max,
                       smoothing_angle, pngname, min=None, max=None):
    """
    Plots 2D half map of signed significance, using Li Ma calculation, masking
    and then smoothing. Saves map as png.

    :param rel_int: String path to fits file with relative intensity
    :param variance: String path to fits file with variance
    :param title: String title of plot
    :param dec_min: float minimum declination for masking
    :param dec_max: float maximum declination for masking
    :param smoothing_angle: float radius for smoothing
    :param pngname: String path to png file that will be saved
    :param min: Optional, float minimum of scale. Default is none.
    :param max: Optional, float minimum of scale. Default is none.

    :return: Saves png file of signed significance.
    """

    # Mask
    masked_combined = mask_map(rel_int, dec_min, dec_max)
    masked_variance = mask_map(variance, dec_min, dec_max)

    # Do np.abs and division on masked array because of hp.UNSEEN
    num = top_hat_smooth(masked_combined, radius=smoothing_angle * DEGREE,
                         average=True)

    masked_np_combined = hp.ma(masked_combined)
    denom = top_hat_smooth(np.abs(masked_np_combined).filled(
        fill_value=hp.UNSEEN), radius=smoothing_angle * DEGREE, average=True)

    masked_np_num = hp.ma(num)
    masked_np_denom = hp.ma(denom)
    # Calculate sign of significance
    neg = masked_np_num / masked_np_denom

    # Calculate Li-Ma Significance
    variance_smooth = top_hat_smooth(masked_variance,
                                     radius=smoothing_angle * DEGREE, average=False)
    masked_np_variance_smooth = hp.ma(variance_smooth)
    significance_llh_iter20 = np.sqrt(np.abs(masked_np_variance_smooth))

    # This makes sure that the significance has the right value and sign
    signed_significance = significance_llh_iter20 * neg
    signed_significance.filled(fill_value=hp.UNSEEN)

    # Plot with specified range and save to directory
    plot_map(signed_significance, title, label='Significance [$\sigma$]',
             half_map=True, dMin=min, dMax=max, filename=pngname)


def tier_one(out_dir, dec_min=-90.0, dec_max=-35.0, smoothing_angle=20,
             range_ri=[0.0, 0.0], range_sig=[0.0, 0.0]):
    """
    Plots relative intensity and significance skymaps for tier one.
    Saves pngs in specified output directory.
    """
    combined = hp.read_map("./fits/combined_t1_iteration20.fits.gz")
    variance = hp.read_map("./fit/significance_t1_iteration20.fits.gz")
    title = "Energy Tier 1: 310 TeV 2011-14"
    if range_ri[0] == 0 and range_ri[1] == 0:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t1_rel_int.png", min=range_ri[0], max=range_ri[1])
    else:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t1_rel_int.png")

    if range_sig[0] == 0 and range_sig[1] == 0:
        li_ma_significance(combined, variance, title, dec_min, dec_max,
                           smoothing_angle, out_dir + "t1_sig.png", min=range_sig[0],
                           max=range_sig[1])
    else:

        li_ma_significance(combined, variance, title, dec_min, dec_max,
                           smoothing_angle, out_dir + "t1_sig.png")


def tier_two(out_dir, dec_min=-90.0, dec_max=-35.0, smoothing_angle=20,
             range_ri=[0.0, 0.0], range_sig=[0.0, 0.0]):
    """
    Plots relative intensity and significance skymaps for tier two.
    Saves pngs in specified output directory.
    """
    combined = hp.read_map("./fits/combined_t2_iteration01.fits.gz")
    variance = hp.read_map("./fits/significance_t2_iteration01.fits.gz")
    title = "Energy Tier 2: 1.1 PeV 2011-14"

    if range_ri[0] == 0 and range_ri[1] == 0:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t2_rel_int.png")
    else:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t2_rel_int.png", min=range_ri[0], max=range_ri[1])

    if range_sig[0] == 0 and range_sig[1] == 0:
        li_ma_significance(combined, variance, title, dec_min, dec_max,
                           smoothing_angle, out_dir + "t2_sig.png")
    else:
        li_ma_significance(combined, variance, title, dec_min, dec_max,
                           smoothing_angle, out_dir + "t2_sig.png", min=range_sig[0],
                           max=range_sig[1])


def tier_three(out_dir, dec_min=-90.0, dec_max=-35.0, smoothing_angle=20,
               range_ri=[0.0, 0.0], range_sig=[0.0, 0.0]):
    """"
    Plots relative intensity and significance skymaps for tier three.
    Saves pngs in specified output directory.
    """
    combined = hp.read_map("./fits/combined_t3_iteration01.fits.gz")
    variance = hp.read_map("./fits/significance_t3_iteration01.fits.gz")
    title = "Energy Tier 3: 2.4 PeV 2011-21"

    if range_ri[0] == 0 and range_ri[1] == 0:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t3_rel_int.png")
    else:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t3_rel_int.png", min=range_ri[0], max=range_ri[1])
    if range_sig[0] == 0 and range_sig[1] == 0:
        li_ma_significance(combined, variance, title, dec_min, dec_max,
                           smoothing_angle, out_dir + "t3_sig.png")
    else:
        li_ma_significance(combined, variance, title, dec_min, dec_max,
                           smoothing_angle, out_dir + "t3_sig.png", min=range_sig[0],
                           max=range_sig[1])


def tier_four(out_dir, dec_min=-90.0, dec_max=-35.0, smoothing_angle=20,
              range_ri=[0.0, 0.0], range_sig=[0.0, 0.0]):
    """
    Plots relative intensity and significance skymaps for tier four.
    Saves pngs in specified output directory.
    """
    combined = hp.read_map("./fits/combined_t4_iteration01.fits.gz")
    variance = hp.read_map("./fits/significance_t4_iteration01.fits.gz")
    title = "Energy Tier 4: 6.6 PeV 2011-21"

    if range_ri[0] == 0 and range_ri[1] == 0:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t4_rel_int.png")
    else:
        relative_intensity(combined, title, dec_min, dec_max, smoothing_angle,
                           out_dir + "t4_rel_int.png", min=range_ri[0], max=range_ri[1])
    if range_sig[0] == 0 and range_sig[1] == 0:
        li_ma_significance(combined, variance, title, dec_min, dec_max,
                           smoothing_angle, out_dir + "t4_sig.png")
    else:
        li_ma_significance(combined, variance, title, dec_min, dec_max,
                       smoothing_angle, out_dir + "t4_sig.png", min=range_sig[0],
                       max=range_sig[1])


def all_tiers(out_dir="./output/", dec_min=-90.0, dec_max=-35.0, smoothing_angle=20):
    """
    Produces 2D relative intensity and significance skymaps for all four tiers
    using optimal scales.

    :param out_dir: String path to output directory. Default is "./output/"
    :param dec_min: float minimum declination for masking. Default is -90.0
    :param dec_max: float maximum declination for masking. Default is -35.0
    :param smoothing_angle: float radius for smoothing in degrees. Default is 20.0.
    
    :return: Saves all 8 skymaps to output directory.
    """
    tier_one(out_dir, dec_min, dec_max, smoothing_angle, range_ri=[-.002, .002], range_sig=[-4.5, 4.5])
    tier_two(out_dir, dec_min, dec_max, smoothing_angle, range_ri=[-.003, .003], range_sig=[-6, 6])
    tier_three(out_dir, dec_min, dec_max, smoothing_angle, range_ri=[-.003, .003], range_sig=[-6, 6])
    tier_four(out_dir, dec_min, dec_max, smoothing_angle, range_ri=[-.005, .005], range_sig=[-4, 4])
    print("DONE")


if __name__ == "__main__":
    # Set up command line options
    parser = argparse.ArgumentParser(description=
                                     "Produce 2D skymaps for specified tiers")
    parser.add_argument("-o", "--output", dest="output", type=str,
                        default="./output/",
                        help="Output directory of .png files")
    parser.add_argument("-d ", "--decmin", dest="decmin", type=float,
                        default=-90.0,
                        help="minimum declination (default = -90.0)")
    parser.add_argument("-D ", "--decmax", dest="decmax", type=float,
                        default=-35.0,
                        help="maximum declination (default = -35.0)")
    parser.add_argument("-s", "--smoothingangle", dest="smooth_angle",
                        type=float, default=20,
                        help="radius of smoothing in degrees")
    parser.add_argument("-a", "--all", action="store_true", dest="all_tiers",
                        default=False, help="Run all tiers")
    parser.add_argument("--tier_one", dest="tier_one", type=float, nargs=2,
                        default=[0, 0], help="Run tier 1. Default none, suggested .002 4.5")
    parser.add_argument("--tier_two", dest="tier_two", type=float, nargs=2,
                        default=[0, 0], help="Run tier 2. Default none, suggested .003 6")
    parser.add_argument("--tier_three", dest="tier_three", type=float, nargs=2,
                        default=[0, 0], help="Run tier 3. Default none, suggested .003 6")
    parser.add_argument("--tier_four", dest="tier_four", type=float, nargs=2,
                        default=[0, 0], help="Run tier 4. Default none, suggested .005 4")
    args = parser.parse_args()

    # If all tiers option, run all tiers
    if args.all_tiers:
        all_tiers(args.output, args.decmin, args.decmax, args.smooth_angle)
    # Else run specified tiers
    else:
        if args.tier_one:
            tier_one(args.output, args.decmin, args.decmax, args.smooth_angle,
                     [-args.tier_one[0], args.tier_one[0]],
                     [-args.tier_one[1], args.tier_one[1]])
        if args.tier_two:
            tier_two(args.output, args.decmin, args.decmax, args.smooth_angle,
                     [-args.tier_two[0], args.tier_two[0]],
                     [-args.tier_two[1], args.tier_two[1]])
        if args.tier_three:
            tier_three(args.output, args.decmin, args.decmax, args.smooth_angle,
                       [-args.tier_three[0], args.tier_three[0]],
                       [-args.tier_three[1], args.tier_three[1]])
        if args.tier_four:
            tier_four(args.output, args.decmin, args.decmax, args.smooth_angle,
                      [-args.tier_four[0], args.tier_four[0]],
                      [-args.tier_four[1], args.tier_four[1]])
