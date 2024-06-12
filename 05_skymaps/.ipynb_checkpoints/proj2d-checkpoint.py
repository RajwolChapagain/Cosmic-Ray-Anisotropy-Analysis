#!/usr/bin/env python

from __future__ import print_function


import numpy as np
import healpy as hp
import os, optparse, re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.container import ErrorbarContainer
from scipy.optimize import curve_fit
from scipy import stats
from mapFunctions import getMap

#rc('xtick', labelsize=18)
#rc('ytick', labelsize=18)

#-------------------------------------------------------



def returnRI(bgmap, datamap, **opts):

    # Setup right-ascension bins
    degree = np.pi / 180
    ramin = opts['ramin'] * degree
    ramax = opts['ramax'] * degree
    rabins = np.linspace(ramin, ramax, opts['nbins']+1)

    # Calculate phi for each pixel
    npix  = len(bgmap)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))

    # Treat masked pixels as zeroes for weighting purposes
    dataweight = np.where((datamap==hp.UNSEEN), 0, datamap)
    bgweight = np.where((bgmap==hp.UNSEEN), 0, bgmap)

    # Bin in right ascension
    data = np.histogram(phi, bins=rabins, weights=dataweight)[0]
    bg   = np.histogram(phi, bins=rabins, weights=bgweight)[0]

    with np.errstate(invalid='ignore'):
        ri = (data - bg) / bg
        sigmay = np.sqrt(data * (bg + 1./20 * data)/bg**3)
    dx = (ramax - ramin)/(2*opts['nbins'])
    ra = np.linspace(ramin+dx, ramax-dx, opts['nbins']) / degree
    sigmax = dx * np.ones(opts['nbins']) / degree

    return (ra, ri, sigmax, sigmay)


def getRelInt(f, **opts):

    # Get background and data maps
    bgmap = getMap([f], mapName='bg', **opts)
    datamap = getMap([f], mapName='data', **opts)
    vals = returnRI(bgmap, datamap, **opts)
    return vals

# Flat line
def flatLine(x, *p):
    return 0*x + p[0]

# Best-fit parameters for flat line
def lineFit(x, y, sigmay):
    popt, pcov = curve_fit(flatLine, x, y, [0], 
                           sigma=sigmay, absolute_sigma=True)
    chi2 = sum((y - flatLine(x, *popt))**2 / sigmay**2)
    return popt, pcov, chi2


# Cosine function with fixed base wavelength of 360 degrees
def multipole(x, *p):
    k = 2*np.pi/360         # Wavenumber (assumes x in degrees)
    l = int(len(p) / 2)    # Multipole number of fit
    return sum([p[2*i] * np.cos((i+1)*k*x - p[2*i+1]) for i in range(l)])

# Best-fit parameters for multipole
def multipoleFit(x, y, l, sigmay):

    # Guess at best fit parameters
    amplitude = (3./np.sqrt(2)) * np.std(y)
    phase     = 0
    p0 = [amplitude, phase] * l

    # Do best fit
    popt, pcov = curve_fit(multipole, x, y, p0, 
                           sigma=sigmay, absolute_sigma=True)
    chi2 = sum((y - multipole(x, *popt))**2 / sigmay**2)

    return popt, pcov, chi2


##=======================================================================##

if __name__ == "__main__":

    # Set up command line options
    usage = "usage: %prog [options] INPUT.fits"
    parser = optparse.OptionParser(usage)
    parser.add_option("-r", "--ramin", dest="ramin", type=float,
            default=0, help="minimum RA")
    parser.add_option("-R", "--ramax", dest="ramax", type=float,
            default=360, help="maximum RA")
    parser.add_option("--rimin", dest="rimin", type=float,
            default=.001, help="minimum relative intensity")
    parser.add_option("--rimax", dest="rimax", type=float,
            default=.001, help="maximum relative intensity")
    parser.add_option("-D", "--decmax", dest="decmax", type=float,
            help="maximum Dec")
    parser.add_option("-d", "--decmin", dest="decmin", type=float,
            help="minimum Dec")
    parser.add_option("-n", "--nbins", dest="nbins", type=int,
            default=24, help="number of bins")
    parser.add_option("-z","--zero", action="store_true", dest="zeroline",
            default=False, help="Draw zero line")
    parser.add_option("-f","--flipra", action="store_true", dest="flipra",
            default=False, help="Flips RA in x axis")
    parser.add_option("-o", "--output", dest="output", default=None,
            help="Output image file name")
    parser.add_option("--multi", dest='multi', type=int,
            default=None, help='Use multipole subtraction')
    parser.add_option('--multiErr', dest='multiErr',
            default=False, action='store_true',
            help='Use amplitude of dipole fit to calculate sys error')
    parser.add_option('--fit', dest='fit',
            type=int,
            help='Show best-fit multipole on plot (1=dipole, 2=quad, etc.)')
    parser.add_option('--offset', dest='offset',
            default=False, action='store_true',
            help='Offset points to avoid overlap')
    parser.add_option('--split', dest='split',
            default=False, action='store_true',
            help='Split input, sharing x-axis')
    parser.add_option("--labels", dest='labels',
            help='Custom label options built-in [configs, method]')
    parser.add_option("-v","--verbose", action="store_true", dest='verbose',
            default=False, help='Optional additional output')
    parser.add_option('--full', dest='full',
            default=False, action='store_true',
            help='Show average behavior of full time range')

    # NOTE: I'd love to clean this up, but I think I just need to add to it. 
    # Need: scale (multiply by RI and show up on axis)
    parser.add_option('-S', '--scale', dest='scale',
            type=int, default=0,
            help='Exponential scale for multiplying y-axis')
    parser.add_option('-L', '--legend', dest='legend',
            default=False, action='store_true',
            help='Display plot legend')
    parser.add_option('--flat', dest='flat',
            default=False, action='store_true',
            help='Show fit parameters for a flat line')
    parser.add_option('-m', '--min', dest='min',
            type=float,
            help='Set plot minimum value')
    parser.add_option('-M', '--max', dest='max',
            type=float,
            help='Set plot maximum value')

    options, args = parser.parse_args()
    opts = vars(options).copy()

    # Default masking behavior
    if not options.decmax:
        opts['mask'] = True

    if options.verbose:
        for key in sorted(opts.keys()):
            print(' --%s %s' % (key, opts[key]))

    # Setup dictionaries for default values
    methods = {'sid':'sidereal', 'solar':'solar', 'anti':'anti-sidereal'}
    methods['ext'] = 'extended-sidereal'
    errDict = {'sid':'anti', 'solar':'ext'}
    p = {}      # Storage for best-fit parameters

    # Plotting setup
    figdims = (8,6)
    if opts['offset']:
        figdims = (17,6)
    axs = None
    if opts['split']:
        fig, axs = plt.subplots(2, sharex=True, sharey=True, figsize=figdims)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    else:
        fig = plt.figure(figsize=figdims)
        ax = fig.add_subplot(111)

    # Setup for referee adjustments
    ramin = opts['ramin']
    ramax = opts['ramax']
    rabins = np.linspace(ramin, ramax, opts['nbins']+1)

    # Systematic error bar setup
    eb = {'edgecolor':None, 'linewidth':0, 'linestyle':None}

    # Y-axis scale
    scale = 10**opts['scale']

    for i, f in enumerate(args):

        # Select plotting axis
        ax = axs[i] if opts['split'] else ax
        
        # Optionally adjust y-axis limits
        if options.min != None and options.max != None:
            ax.set_ylim(options.min, options.max)

        # Labeling
        basename = os.path.basename(f)[:-5]
        label = basename
        # Specialty labeling based on configuration or time frame
        if options.labels == 'configs':
            label = re.findall('IC86-\d{4}', f)[-1]
        # Find time frame, stripping leading/trailing characters
        #method = re.findall('|'.join([f'_{k}[_|\.]' for k in methods]), f)[-1]
        method = "sidereal"
        if options.labels == 'method':
            label = method

        # Possible file for systematic uncertainties
        errfile = ''
        if method in errDict.keys():
            errfile = f.replace(method, errDict[method])

        # Get RA, relint, and errors
        ra, ri, sigmax, sigmay = getRelInt(f, **opts)
        ri *= scale
        sigmay *= scale
        if opts['offset']:
            npoints = len(args)
            # Offset for time-dependent plot
            dx = (rabins[1] - rabins[0]) / float(npoints+1)
            ra = rabins[:-1] + (npoints-i)*dx

        # Plot 1-D relative intensity projection
        l = ax.errorbar(ra, ri, xerr=0*sigmax, yerr=sigmay,
                    marker='.', fmt='.',
                    capsize=4, label=label, linewidth=2, markersize=8, mew=0)
        if options.verbose:
            print("sigma: ", basename, np.std(ri))
            print("max: ", basename, np.max(np.abs(ri)))

        # Optional best fit(s)
        if options.fit != None:

            # Calculate best-fit line and uncertainties
            popt, pcov, chi2 = multipoleFit(ra, ri, options.fit, sigmay)
            eopt = np.sqrt(np.diag(pcov))

            # Plot a smoothed version of the applied fit
            tOpt = {'color':'blue'}
            smooth_x = np.linspace(ra[0], ra[-1], ra.size*20)
            ax.plot(smooth_x, multipole(smooth_x, *popt), **tOpt)

            # Calculate degrees of freedom and p-value
            ndof = ra.size - popt.size
            pvalue = 1 - stats.chi2.cdf(chi2, ndof)

            # Text locations and values
            ymax = ax.get_ylim()[-1]
            delta = ymax/10
            x0 = 100    # x-location (in degrees)
            y0 = 0 if method=='sid' else ymax
            # Adjust best-fit amplitude and phase to positive values
            amp = popt[0]
            phase = popt[1] * 180/np.pi
            if amp < 0:
                phase -= 180
                amp *= -1
            while phase < 0:
                phase += 360
            phase_err = eopt[1] * 180/np.pi

            # Write to plot
            tOpt.update({'fontweight':'bold'})
            info = {}
            info[0] = f'A = {amp:0.2f} $\pm$ {eopt[0]:0.2f}'
            info[1] = f'$\phi$ = {phase:0.1f} $\pm$ {phase_err:0.1f}'
            info[2] = f'$\chi^2$/NDOF = {chi2:0.0f} $/$ {ndof:0.0f}'
            if method != 'sid':
                info[3] = f'Probability = {pvalue:0.4f}'
            for i in info.keys():
                ax.text(x0, y0-(i+1)*delta, info[i], **tOpt)

        if options.flat:

            # Calculate and plot flat line
            popt, pcov, chi2 = lineFit(ra, ri, sigmay)
            eopt = np.sqrt(np.diag(pcov))
            tOpt = {'color':'orange'}
            ax.plot(ra, flatLine(ra, *popt), **tOpt)

            # Calculate degrees of freedom and p-value
            ndof = ra.size - popt.size
            pvalue = 1 - stats.chi2.cdf(chi2, ndof)

            # Text locations and values
            ymax = ax.get_ylim()[-1]
            delta = ymax/10
            x0 = 350    # x-location (in degrees)
            y0 = 0 if method=='sid' else ymax

            # Write to plot
            tOpt.update({'fontweight':'bold'})
            info = {}
            info[0] = f'$\chi^2$/NDOF = {chi2:0.0f} $/$ {ndof:0.0f}'
            info[1] = f'Probability = {pvalue:0.4f}'
            for i in info.keys():
                ax.text(x0, y0-(i+1)*delta, info[i], **tOpt)

        # Additional plotting options
        if not os.path.isfile(errfile) and errfile!='':
            print(f'Error file {errfile} not found! Cannot calculate syserr')

        if os.path.isfile(errfile):

            ra_err, ri_err, sigx, sigy = getRelInt(errfile, **opts)
            ri_err *= scale
            sigy *= scale

            # Ideal scenario - fit cosine function and take amplitude
            # Not used because anti- and extended-sidereal not sinusoidal
            if options.multiErr:
                popt, pcov, chi2  = multipoleFit(ra_err, ri_err, 1, sigy)
                syserr = popt[0]
            else:
                #syserr = np.abs(ri_err).max()  # Six-year method
                syserr = np.sqrt(np.mean(ri_err**2))

            box = dx if opts['offset'] else 10
            patches = [mpl.patches.Rectangle([ra[j]-box/2, ri[j]-syserr], box,
                    2*syserr, **eb) for j in range(len(ra))]
            cln = PatchCollection(patches, cmap=mpl.cm.jet,
                    alpha=0.5, facecolor=l[0].get_color())
            ax.add_collection(cln)

        # Add dashed lines to separate right ascension bins
        if opts['offset']:
            for raval in rabins[1:-1]:
                ax.axvline(x=raval, color='k', ls=':')

        # Axes labels
        ax.set_xlabel(r"Right Ascension $[^{\circ}]$",fontsize=14)
        if opts['split']:
            ylabel = r'$\Delta N/\langle N \rangle$'
            x0, y0 = -0.04, 0.5
            if opts['scale'] != 0:
                ylabel += fr'$\; (\times 10^{{{-opts["scale"]}}})$'
                x0, y0 = 0.05, 0.5
            fig.text(x0, y0, ylabel, fontsize=14,
                    va='center', rotation='vertical')
        else:
            ylabel = r'$\Delta N/\langle N \rangle$'
            if opts['scale'] != 0:
                ylabel += fr'$\;(\times 10^{{{-opts["scale"]}}})$'
            ax.set_ylabel(ylabel,fontsize=14)
        #ax.grid()

    # Optionally include average of whole dataset
    axs = axs if opts['split'] else [ax]
    for ax in axs:
        if options.full:
            totalMap = re.sub('IC86-\d{4}_','IC86_', f)
            ra, ri, sigmax, sigmay = getRelInt(totalMap, **opts)
            ri *= scale
            sigmay *= scale
            totalErr = re.sub('IC86-\d{4}_','IC86_', errfile)
            if os.path.isfile(totalErr):
                print(f'Error file {os.path.basename(errfile)} found.')
                ra_err, ri_err, sigx, sigy = getRelInt(totalErr, **opts)
                ri_err *= scale
                sigy *= scale
                if options.multiErr:
                    popt, pcov, chi2 = multipoleFit(ra_err, ri_err, 1, sigy)
                    syserr = popt[0]
                else:
                    #syserr = np.abs(ri_err).max()  # Six-year method
                    syserr = np.sqrt(np.mean(ri_err**2))
                box = dx * len(args)
                patches = [mpl.patches.Rectangle([ra[j]-box/2, ri[j]-syserr],
                        box, 2*syserr, **eb) for j in range(len(ra))]
                cln = PatchCollection(patches, alpha=0.4, facecolor='gray')
                ax.add_collection(cln)

        # Show zero line
        if options.zeroline:
            xzero = np.arange(0, 360, 1)
            yzero = 0 * xzero
            ax.plot(xzero,yzero,linewidth=1.5,linestyle='--',color='black')

        # Adjust x-axis limits
        ax.set_xlim(options.ramax, options.ramin)
        if options.flipra:
            ax.set_xlim(options.ramin, options.ramax)

        # Show legend
        if opts['legend']:
            leg = ax.legend(loc='lower right')

    plt.draw()
    plt.savefig(options.output, dpi=300, bbox_inches='tight')
    

