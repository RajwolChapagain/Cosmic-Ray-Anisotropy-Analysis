import os, fileinput, sys


def runSim(particle, fnum):
	for line in fileinput.input('steeringcard_sim',inplace=1):
		arg_string = particle + ' ' + str(fnum)
		if 'arguments' in line:
			line = line.replace(line, 'arguments = '+arg_string+'\n')
		sys.stdout.write(line)
	os.system('condor_submit steeringcard_dstsim')
	print('running for', particle, fnum)


#p 12360
#He 12630 Me
#O 12631 Me
#Fe 362

runSim("p", 12360)
runSim("He", 12630)
runSim("O", 12631)
runSim("Fe", 12362)

