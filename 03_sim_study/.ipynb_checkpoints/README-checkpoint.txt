Instructions to produce plots based on 2012 Monte Carlo simulation of IceTop

Step 1: Processing sim hdf5 files into .npy arrays
* In icetray, run 'python 001_nstation_binning.py /data/user/@USER_DIR@' with
  the user directory where the simulation hdf5 files are located as an argument.
  (Ignore the warnings)

 - This will read the hdf5 files, bin them by the number of stations triggered,
   and save the data in accessible numpy arrays in ./arrays. 

 - There are three options. Defaults are applied for our analysis.
   * --it73: Integer from 0 to 2, optional
        Index of IT73 Quality Cuts Permutation
        0 -> No quality cuts [default, what we used in our analysis]
        1 -> IceTop_reco_succeeded is applied
        2 -> All IT73 Quality Cuts are applied
   * -o, --output: String, optional 
        Path to output directory where the .npy arrays will be saved
        Default is "./arrays". If changed, edit line 3 in cell 2 of 
        002_loading_graphs.ipynb. 
   * -e, --energies: String of lists (default: "[[3, 5], [5, 9], [9, 16]]")
        Each element is a list of length two, containing the low
        and high bounds of each energy bin. E.g. energies = [[3,5]] means a
        low energy bin for 3 <= number of stations < 5, and a high energy bin
        for 5 <= number of stations

  
Step 2: Producing plots.
* To produce the plots, run all cells in the jupyter notebook 
  '002_Binning_Graphs.ipynb' 
 - Kernel: py3-v4.1.1.:combo/V01-01-06

 - These plots are created using functions in ang_res_funcs.py. The font size,
   colour, etc. can be modified here. 