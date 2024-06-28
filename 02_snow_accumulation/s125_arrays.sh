#!/bin/bash

#eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
#eval '/cvmfs/icecube.opensciencegrid.org/standard/icetray-start'
#eval cvmfs
#eval combo


TMPDIR=$(mktemp -d)
cp s125_arrays.py $TMPDIR

cd $TMPDIR
python s125_arrays.py $1

