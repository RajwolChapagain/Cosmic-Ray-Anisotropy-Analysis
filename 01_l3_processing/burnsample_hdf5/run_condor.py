import os, fileinput, sys

##run this code in condor shell to submit condor job for each month of data (run through i3 convert py file)
#this can be used for creating burnsample hdf5 and root files
for run_year in range(2011, 2022):
    for year in [run_year]:
        for month in range(4,13):
            for line in fileinput.input('steeringcard_dstburnsample',inplace=1):
                arg_string = str(run_year)+' '+str(year)+' '+str(month)
                if 'arguments' in line:
                    line = line.replace(line,'arguments = '+arg_string+'\n')
                sys.stdout.write(line)
	    os.system('condor_submit steeringcard_dstburnsample') 
	    print('running for', run_year, year, month)


    for year in [run_year+1]:
        for month in range(1, 8):
            for line in fileinput.input('steeringcard_dstburnsample',inplace=1):
	        arg_string = str(run_year)+' '+str(year)+' '+str(month)
		if 'arguments' in line:
		    line = line.replace(line,'arguments = '+arg_string+'\n')
		sys.stdout.write(line)
	    os.system('condor_submit steeringcard_dstburnsample')
	    print('running for', run_year, year, month)

