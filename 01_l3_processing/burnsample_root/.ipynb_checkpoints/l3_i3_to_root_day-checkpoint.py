from icecube import icetray, dataclasses, simclasses, dataio, tableio, toprec, astro

from icecube.rootwriter import I3ROOTWriter
from I3Tray import I3Tray
from icecube.icetop_Level3_scripts.functions import count_stations

import sys 
import numpy as np
from argparse import ArgumentParser
from glob import glob

def add_IceTop_quality_cuts(frame):
        passed = all(frame['IT73AnalysisIceTopQualityCuts'].values())
        #print "test " + str(passed)
        frame['passed_IceTopQualityCuts'] = icetray.I3Bool(passed)


def add_nstations(frame, pulses='IceTopHLCSeedRTPulses_SnowUnAttenuated'):
    nstation = count_stations(dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulses))
    frame['NStations'] = icetray.I3Int(nstation)

#def add_get_equatorial(frame, 


parser = ArgumentParser(description=__doc__)
parser.add_argument('infile', nargs='*')
#parser.add_argument('outfile')
opts = parser.parse_args()

outfile_path = '' #update to correct outfile path

run_year = int(sys.argv[1])
year = int(sys.argv[2])
month = int(sys.argv[3])
day = int(sys.argv[4])


physics = glob("/data/ana/CosmicRay/IceTop_level3/exp/IC86."+str(run_year)+"_pass2_v*/"+str(year)+"/"+str(month).zfill(2)+str(day).zfill(2)+"/*/Level3_IC86."+str(run_year)+"*_Subrun*.i3.*")
#glob() the list of files from the disk


physics.sort()


subevent = 'IceTopSplit'

temp_outfile = 'l3_data_run_config_'+str(run_year)+'_year_'+str(year)+'_'+str(month).zfill(2)+str(day).zfill(2)+'.root'
if len(physics)!=0:
	tray = I3Tray()
	#tray.AddModule('I3Reader', 'reader', filenamelist=opts.infile)
	tray.AddModule('I3Reader', 'reader', DropBuffers= False, filenamelist=physics)

	tray.Add(add_nstations,
		 If=lambda frame: 'IceTopHLCSeedRTPulses_SnowUnAttenuated' in frame)

	tray.Add(add_IceTop_quality_cuts,
		 If=lambda frame: 'IT73AnalysisIceTopQualityCuts' in frame)

##	tray.Add(add_get_equatorial,

	tray.AddSegment(I3ROOTWriter, 'scribe',
	    #Output=opts.outfile,
	    Output=temp_outfile,
	    Keys=[
		  'I3EventHeader',
		  'QFilterMask',
		  'Laputop',
		  'LaputopParams',
		  'LaputopSmall',
		  'LaputopSmallParams',
		  'ShowerPlane',
		  'ShowerPlaneParams',
		  'IT73AnalysisIceTopQualityCuts',
		  'NStations'
		  ],

	    Types=[],
	    SubEventStreams=[subevent],
	    #DropOrphanStreams=[icetray.I3Frame.DAQ],
	    #streams=[icetray.I3Frame.Physics],
	)


	tray.Execute()																															

