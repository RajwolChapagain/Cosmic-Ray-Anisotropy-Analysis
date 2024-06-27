import os, fileinput, sys

##run this code in condor shell to submit condor job for each month of data (run through i3 convert py file)
#this can be used for creating burnsample hdf5 and root files
for run_year in range(2011, 2022):
    for year in [run_year]:
        for month in range(4,13):
            for day in range(1, 32):
               	for line in fileinput.input('steeringcard_dstroot',inplace=1):
                    arg_string = str(run_year)+' '+str(year)+' '+str(month)+' '+str(day)
                    if 'arguments' in line:
                        line = line.replace(line,'arguments = '+arg_string+'\n')
                    sys.stdout.write(line)
		os.system('condor_submit steeringcard_dstroot')
		print('running for', run_year, year, month, day)


    for year in [run_year+1]:
        for month in range(1, 6):
	    for day in range(1, 32):
	        for line in fileinput.input('steeringcard_dst',inplace=1):
		    arg_string = str(run_year)+' '+str(year)+' '+str(month)+' '+str(day)
		    if 'arguments' in line:
		        line = line.replace(line,'arguments = '+arg_string+'\n')
		    sys.stdout.write(line)
		os.system('condor_submit steeringcard_dstroot') 
		print('running for', run_year, year, month, day)

