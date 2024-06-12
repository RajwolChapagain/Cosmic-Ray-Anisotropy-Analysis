#!/bin/bash

input='./fits/'
output='./output/'
fit='--gaussian'
#range='-M 12 -m -12'
range=''
python 1d_proj.py $input"combined_t1_iteration20.fits.gz" -o $output"1d_t1_gauss.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 1: 310 TeV 2011-14"
python 1d_proj.py $input"combined_t2_iteration20.fits.gz" -o $output"1d_t2_gauss.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 2: 1.1 PeV 2011-14"
python 1d_proj.py $input"combined_t3_iteration20.fits.gz" -o $output"1d_t3_gauss.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 3: 2.4 PeV 2011-21"
python 1d_proj.py $input"combined_t4_iteration20.fits.gz" -o $output"1d_t4_gauss.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 4: 6.6 PeV 2011-21"

