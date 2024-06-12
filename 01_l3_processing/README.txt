Instructions to process level3 files into .hdf5 and .root files for later analysis

Step 0: Editing the code
- Edit the output directories in the condor wrapper shell scripts in each directory.
- Edit the output, log, and scratch destination directories in the steeringcards in each directory.

Step 1: Running
- We produce .hdf5 files of the burnsample in order to do the snow accumulation study. 
 * Run ./burnsample_hdf5/run_condor.py in icetray in condor.
- We produce .hdf5 files of the simulation in order to produce simulation plots.
 * Run ./sim_hdf5/run_stuff.py in icetray in condor.
- We produce .root files of the burnsample in order to do the anisotropy analysis.
 * Run ./burnsample_root/run_condor.py in icetray in condor.