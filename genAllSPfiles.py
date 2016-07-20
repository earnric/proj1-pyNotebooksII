#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 24 Feb 2016
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and outputs particle data
#
# Revision history
#

# ##########################################################
# Output the star particle data for all stars in the snap
# ##########################################################
def outputStarData(z,s):
    print ("Saving data...")
    np.savetxt("spLoc_%05.2lf.txt"%z,s.s['pos'],header="x\ty\tz",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spMass_%05.2lf.txt"%z,s.s['mass'],header="mass",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spZ_%05.2lf.txt"%z,s.s['zsolar'],header="Z_solar",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spPPF_%05.2lf.txt"%z,s.s['ppf'],header="ppf",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spPZ_%05.2lf.txt"%z,s.s['pzsolar'],header="pz_solar",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spBT_%05.2lf.txt"%z,s.s['tform'].in_units('yr'),header="birth_time_yr",comments='') # comments='' gets rid of "#" in output
    print ("Star Particle Data Files saved for %05.1lf"%z)
    return


# ##########################################################
# ##########################################################
##
## Main program
##
# ##########################################################
# ##########################################################
#%matplotlib inline
import os
import re # Regular expression matching
import numpy as np
import pynbody
import sys, os, glob, pickle, gc
import math
#import mmap

pynbody.ramses.multiprocess_num = 8
pynbody.config['number_of_threads'] = 32

# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
files = [f for f in os.listdir('.') if re.match(r'output_000[0-9][0-9]', f)]
files =[
    # "output_00007",
    # "output_00008",
    # "output_00011",
    "output_00016",
    # "output_00020",
    # "output_00026",
    # "output_00033",
    # "output_00037",
    # "output_00043",
    # "output_00050",
    # "output_00058",
    # "output_00066",
    # "output_00068",
    "output_00073",
    # "output_00084",
    # "output_00097",
    # "output_00107",
    "output_00121",
    # "output_00136",
    # "output_00152",
    # "output_00169",
    "output_00191",
    "output_00214",
    "output_00237"
    ]
files.sort()
print ("Files to process %d"%len(files))

np.set_printoptions(precision=3,suppress=True)

# Loop and load data
# If we wanted to just loop over filename - "for file in files:"
for indx, file in enumerate(files):
    gc.collect()
    print ("Processing file %s" % file)
    s = pynbody.load("./"+file)
    s['pos'] -= 0.5
    s.physical_units();
    
    z = 1/s.properties['a']-1
    print ("Redshift = %.2lf" % z)
    boxsizestring = "%.2lf kpc" % s.properties['boxsize'].in_units('kpc')
    print ("Boxsize @ this redshift %s"%boxsizestring)
    print ("aexp %.3lf" % s.properties['a'])
    bs = float(s.properties['boxsize'])

    z_crit = 2.0e-7 # Mass fraction
    print ("Normalize Z data...")
    ## s.g['metal'][s.g['metal']<1e-10]    = 1e-10 # Since we divide by Z below, don't use 0.0
    ## s.g['pgf'][s.g['pgf']>(1.0-1e-6)]  = 1.0
    ## s.g['pzf'][s.g['pzf']<1e-10]        = 1e-10
    ## s.g['zsolar'] = s.g['metal'] * 50.0         # Solar units
    ## s.g['pzsolar'] = s.g['pzf'] * 50.0          # Solar units

    s.s['metal'][s.s['metal']<1e-10]    = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-6)]   = 1.0
    s.s['pzf'][s.s['pzf']<1e-10]        = 1e-10
    s.s['zsolar'] = s.s['metal'] * 50.0         # Solar units
    s.s['pzsolar'] = s.s['pzf'] * 50.0          # Solar units

    # Fix birth times
    print (s.derivable_keys()[1:5])
    pynbody.analysis.ramses_util.get_tform(s,"/Users/earnric/bin/part2birth");

    # ##########################################################
    # Output SP data
    outputStarData(z,s)
 
    del s
    gc.collect()
    gc.collect()

print ("Done")
