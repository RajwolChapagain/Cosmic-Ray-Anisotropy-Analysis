#!/bin/bash

#eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
#eval '/cvmfs/icecube.opensciencegrid.org/standard/icetray-start'
#eval cvmfs
#eval combo


TMPDIR=$(mktemp -d)
cp l3_i3_to_root_day.py $TMPDIR

cd $TMPDIR
python l3_i3_to_root_day.py $1 $2 $3 $4

mv l3_data* /data/user/@USER_DIR@/burnsamplemaps


