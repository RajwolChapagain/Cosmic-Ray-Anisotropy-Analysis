Instructions to produce 2D Skymaps and 1D Projections of Relative Intensity and Significance.

Step 0: Setting up the input directory.
* The 2Dskymaps python script assumes that the fits files follow the below
  schema:
 - The CR_IceTop_64_360_iteration20.fits.gz files produced in the previous
   step are copied to the fits directory here and saved as 
   "combined_tX_iteration20.fits.gz" where X is the numeral for each tier
   (e.g. combined_t1_iteration20.fits.gz)

 - The significance_IceTop_64_360_iteration20.fits.gz files produced in the
   previous step are copied to the fits directory here and saved as
   "significance_tX_iteration20.fits.gz" where X is the numeral for each tier
   (e.g. significance_t1_iteration20.fits.gz)

* You should have one combined fits file and one significance fits file for
  each tier, 8 files total in ./fits/

Running the code
Step 1: Producing 2D Skymaps
* In icetray, run "python 2dskymaps.py [opts]"
 - You MUST specify which tiers to produce skymaps for. Tiers can be listed
   (e.g. --tier_one --tier_three) or you can use "-a" to run all tiers.
  - IF you run specific tiers, you can specify the ranges for the relative
    intensity and significance maps by listing two floats after the option. 
    0 will result in the map automatically adjusting the bounds based on the data.
    The suggested bounds, which are automatically applied with "-a" are in the help.

 - The output directory default is ./output/. This can be changed.

 - The decmin, decmax and smoothing angle defaults are what were used for the
   burnsample analysis (-90.0, -35.0 and 20 respectively). These can be changed
   with "-d", "-D" and "-s" respectively.
 - 
Step 2: Producing 1D Projections of Relative Intensity vs Right Ascension
* In cobalt, with CVMFS py3-v4.2.1., run one of the following:
 - ./alltiers.sh
  * This runs the script for all four tiers. The input and output directories
    are the defaults: ./fits/ and ./output/. The fit can be specified:
    --gaussian for a gaussian fit, --flat for a flat fit, --fit 1 and --fit 2
    for first- and second- order harmonics. A range can also be specified
    for all four tiers for uniformity. 
 - python 1d_proj.py INPUT.fits [opts]
  * This script produces a 1D projection of relative intensity based on a fits
    file and per the options. The command will look like:
    python 1d_proj.py /data/user/@USER_DIR@/combined_t1_iteration20.fits.gz -o ./output/1d_t1_gauss.png --gaussian -M 2 -m -2 -d -90 -D -35 -z -S 3 --title "Tier 1: 310 TeV 2011-14"
