#!/bin/bash

#eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
#eval '/cvmfs/icecube.opensciencegrid.org/standard/icetray-start'
#eval cvmfs
#eval combo


TMPDIR=$(mktemp -d)
cp l3_i3_burnsample_hdf.py $TMPDIR

cd $TMPDIR
python l3_i3_burnsample_hdf.py $1 $2 $3

mv l3_data* /data/user/@USER_DIR/burnsamplehdf5


