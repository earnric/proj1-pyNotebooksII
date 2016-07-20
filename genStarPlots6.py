#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 09 Sep 2015
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and computes the various
#  pynbody as well as star-location plots.
#
# Method:
#  Sample 3 sps from the s.s collection and plot Z, PGF,
#  Temp, Rho and locations for each... 
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
# Output the data needed to generate the Mathematica histograms
# ##########################################################
def outputHistogramData(z,i,snBlast):
    bins    = np.logspace(-10, 0, 51) # Log bins for histogram data 
    psm     = np.zeros(len(bins)-1) # Pristine Stellar mass in that bin
    tsm     = np.zeros(len(bins)-1) # total mass
    primsm  = np.zeros(len(bins)-1) # Primordial stellar mass

    # For histograms, if pzf close to 0, make it 0
    # If we don't do this, we end up with pzf/Z = 1.0 ... so it would look like
    # everything was pristine!
    temp = s.s
    temp['pzf'][temp['pzf'] < 1e-10] = 0.0
    # snBlast.s['pzf'][snBlast.s['pzf'] <= 1e-10] = 0.0

    # By using zsolar, we are binning in units of solar Z.
    for indx2,j in enumerate(bins):
        if indx2 < len(bins)-1:
            cond = (snBlast.s['zsolar'] >= j) & (snBlast.s['zsolar'] < bins[indx2+1]) # This selects for sp's in the bin
            psm[indx2] = np.sum(snBlast.s['ppf'][cond] * snBlast.s['mass'][cond])
            tsm[indx2] = np.sum(snBlast.s['mass'][cond])
            # For sp's that are in our bin (above):
            #   Compute the polluted fraction (1-ppf)
            #   Compute the fraction of pristine metals: pzf/Z
            #   Compute the mass of polluted stars that are polluted only by pristine metals
            primsm[indx2] = np.sum((1.0-snBlast.s['ppf'][cond]) * 
                                   (temp['pzf'][cond] / snBlast.s['metal'][cond]) * snBlast.s['mass'][cond])
    # Save the total masses
    np.savetxt("primordSM_z%.1lf-%i.txt"%(z,i),primsm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("pristSM_z%.1lf-%i.txt"%(z,i),psm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("totSM_z%.1lf-%i.txt"%(z,i),tsm,comments='') # comments='' gets rid of "#" in output
    del temp
    return

# ##########################################################
# Generate and save images... 
# ##########################################################
def outputImages(z,i,snBlast):
    # Generate sph images for out box...
    with pynbody.transformation.translate(snBlast,coords):
        try:
            fileOut="img_Z-z=%.1lf-%i.png"% (z,i)
            titleStr = "$Z_{\odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="zsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,
                      log=True, vmax=1.0, vmin=1e-7, qtytitle="$log_{10}\;Z_{\odot}$", approximate_fast=False,title=titleStr,
                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            fileOut="img_PGF-z=%.1lf-%i.png"% (z,i)
            titleStr = "PGF - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="pgf",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,
                      log=True, vmax = 1.0, vmin=1e-7, qtytitle="$log_{10}\; PGF$",approximate_fast=False,title=titleStr,
                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            fileOut="img_PZ-z=%.1lf-%i.png"% (z,i)
            titleStr = "$Z_{P, \odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="pzsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,
                      log=True, vmax=1.0, vmin=1e-7, qtytitle="$log_{10}\; Z_{P, \odot}$",approximate_fast=False,title=titleStr,
                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            fileOut="img_Density-z=%.1lf-%i.png"% (z,i)
            titleStr = "Density - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="rho",width=smallbox,cmap="terrain", denoise=True ,av_z=False,
                      log=True, approximate_fast=False,title=titleStr,
                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            fileOut="img_Temp-z=%.1lf-%i.png"% (z,i)
            titleStr = "Temp - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="temp",width=smallbox,cmap="RdYlBu", denoise=True ,av_z=False,
                      log=True, approximate_fast=False,title=titleStr,
                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            fileOut="img_vt-z=%.1lf-%i.png"% (z,i)
            titleStr = "$v_{t}$ - z = %.1lf" % z # + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="tv",width=smallbox,cmap="RdYlBu", denoise=True ,av_z=False,
                      log=True, qtytitle="$log_{10}\; v_{t}$",approximate_fast=False,title=titleStr,
                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            gc.collect()
        except:
            print ("EXCEPTION while processing file: %s"%file)
            print ("Plot title: %s"%titleStr)
            gc.collect()
            pass
    return

# ##########################################################
# Generate and save stellar mass totals for z
# ##########################################################
def outputMassTots(file,s,data,indx):
    # Make a local copy of the stars, set pzf to 0.0
    # for small values.
    temp = s.s
    temp['pzf'][temp['pzf'] < 1e-10] = 0.0
    #s.s['pzf'][s.s['pzf'] <= 1e-10] = 0.0

    # Get total stellar mass
    totalStarMass = np.sum(s.s['mass'].in_units('Msol'))
    print ("Total stellar mass is %.2e Msol @ z~%.2lf" % (totalStarMass,z))
    
    # Total POP3 star mass
    totalPop3StarMass = np.sum(s.s['mass'].in_units('Msol') * s.s['ppf'])
    print ("Total POPIII stellar mass is %.2e Msol @ z~%.2lf" % (totalPop3StarMass,z))

    # Total of only totally pristine star particles ***** Entire SP must be < z_crit!! ****** 
    # Note that we have to use mass-fraction "metal" here... 
    totalPop3SubcritZStarMass = np.sum(s.s['mass'][s.s['metal'] < z_crit])
    print ("Total sp mass for spZ < Z_crit is %.2e Msol @ z~%.2lf" % (totalPop3SubcritZStarMass,z))
    
    # Total polluted stars
    totalPolStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0-s.s['ppf']))
    print ("Total polluted stellar mass is %.2e Msol @ z~%.2lf" % (totalPolStarMass,z))

    # Compute the mass of stars' primordial metals
    # 1 - ppf is the polluted fraction of stars, by mass
    # pzf/Z is then the fraction of primordial metals
    # (1-ppf) * pzf / Z is hence the fraction of stars polluted only by primordial metals
    totalPriStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0 - s.s['ppf']) * temp['pzf']/ s.s['metal'])
    print ("Total primordial stellar mass is %.2e Msol @ z~%.2lf" % (totalPriStarMass,z))

    totalGasMass  = np.sum(s.g['mass'].in_units('Msol'))
    print ("Total gas mass is %.2e Msol @ z~%.2lf" % (totalGasMass,z))

    totalPristGasMass = np.sum(s.g['mass'].in_units('Msol') * s.g['pgf'])
    print ("Total pristing gas mass is %.2e Msol @ z~%.2lf" % (totalPristGasMass,z))

    # Compute total star mass for stars NOT made of primordial Z
    totNonPrimordStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0 - s.s['ppf']) * (1.0 - temp['pzf']/s.s['metal']))
    print ("Total non-primordial stellar mass is %.2e Msol @ z~%.2lf" % (totNonPrimordStarMass, z))

    data[indx] = [z, min(s.s['tform'].in_units('yr')), max(s.s['tform'].in_units('yr')),
                  totalStarMass,totalPop3StarMass,totalPolStarMass,totalPriStarMass,
                  totalGasMass,totalPristGasMass,totalPop3SubcritZStarMass,totNonPrimordStarMass]
    
    np.savetxt("%s_new-data3Mpc.txt"%file,data[0:indx],comments='',header=headerStr) # comments='' gets rid of "#" in output
    del temp
    return data

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
import matplotlib as mpl
mpl.use('Agg') # Needed if you don't have an X-server
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import pynbody.plot as pp
import pynbody.plot.sph as sph
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import sys, os, glob, pickle, gc
import math
#import mmap

mpl.rcParams['figure.figsize'] = (12,8)
mpl.rcParams['font.size'] = 18
pynbody.ramses.multiprocess_num = 32
pynbody.config['number_of_threads'] = 128

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
    # "output_00073",
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

key = np.empty([len(files),2]) # a map of z to boxsize, 

data      = np.empty([len(files),11])
headerStr = 'z\ttStart\ttEnd\ttotStarMass\ttotPop3StarMass\ttotPollStarMass\ttotPrimordStarMass\ttotGasMass\ttotPristGasMass\ttotSubcritStarMass\ttotNonPrimordStarMass'

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
    s.g['metal'][s.g['metal']<1e-10]    = 1e-10 # Since we divide by Z below, don't use 0.0
    s.g['pgf'][s.g['pgf']>(1.0-1e-6)]  = 1.0
    s.g['pzf'][s.g['pzf']<1e-10]        = 1e-10
    s.g['zsolar'] = s.g['metal'] * 50.0         # Solar units
    s.g['pzsolar'] = s.g['pzf'] * 50.0          # Solar units

    s.s['metal'][s.s['metal']<1e-10]    = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-6)]  = 1.0
    s.s['pzf'][s.s['pzf']<1e-10]        = 1e-10
    s.s['zsolar'] = s.s['metal'] * 50.0         # Solar units
    s.s['pzsolar'] = s.s['pzf'] * 50.0          # Solar units

    # Fix birth times
    print (s.derivable_keys()[1:5])
    pynbody.analysis.ramses_util.get_tform(s,"/home1/02744/earnric/bin/part2birth");

    # ##########################################################
    # Output SP data
    outputStarData(z,s)
    # Save a list of redshift and the boxsize for the files we're processing
    key[indx] = ["%.2lf"%z,"%.2lf"%s.properties['boxsize'].in_units('kpc')]

    # ##########################################################
    # Output mass totals for Mathematica overall histogram
    data = outputMassTots(file,s,data,indx)

    sbox = 40.0 / (1.0 + z) * 0.71 # 40 kpc comoving box
    print ("40 kpc comoving Small box @ z=%.2lf is %.2lf"%(z,sbox))

    step = len(s.s)/3 # Note integer division
    if step == 0: step = 1
    for i in range(0,len(s.s),step):
        x,y,zz = - s.s['pos'][i] # get a star
        print ("Star offset [%.2lf %.2lf %.2lf]"%(x,y,zz))
        print ("Boxsize @ this z ", bs)
        # Ensure we center such that we can depict a 40 kpc box
        if (abs(x) > bs/2.0 - sbox/2.0):
            x = np.sign(x) * (bs/2.0 - sbox/2.0) # closest x coord to sp
        if (abs(y) > bs/2.0 - sbox/2.0):
            y = np.sign(y) * (bs/2.0 - sbox/2.0) # closest y coord to sp
#        if (abs(zz) > bs/2.0 - sbox/2.0): # Commented out so we get the exact z-coord of the sp
#            zz = np.sign(zz) * (bs/2.0 - sbox/2.0) # we could go off the end... ok
                    
        # Create the size of our box for sph image plots
        # May not be centered on sp - we have adjusted x,y such that we
        # don't go 'off the edge' for a 'smallbox' sized image
        coords= [x,y,zz] # NOTE WE ARE TRANSLATING... coords is -1 * location of sp of interest
        if sbox > 1.0: smallbox= str(sbox) + ' kpc'
        else: smallbox= str(sbox*1000.0) + ' pc'
        print ("Plotsize: ",smallbox)

        # ##########################################################
        # Save the index and coordinates of the center of our image...
        np.savetxt("z%05.2lf_SpCoord_%i.txt"%(z,i), -1.0*np.array(coords),comments='') # comments='' gets rid of "#" in output
        # ##########################################################

        rx,ry,rz = -1 * np.array([x,y,zz]) # For filters, we need the actual coords, not translation coords.
        print ("Cuboid center [%.2lf %.2lf %.2lf]"%(rx,ry,rz))
        # Note that the thickness in z is only 1/2 that of the other dimensions
        snBlast = s[pynbody.filt.Cuboid(str((rx-sbox/2.0)) + " kpc", str((ry-sbox/2.0)) + " kpc",str((rz-sbox/4.0)) + " kpc",
                                        str((rx+sbox/2.0)) + " kpc", str((ry+sbox/2.0)) + " kpc",str((rz+sbox/4.0)) + " kpc")]
        
        # Save a DM image of the entire box at this z.
        ## try:
        ##     titleStr = "DM - z = %.1lf"%z
        ##     fileOut  = "DM-z=%.lf.png"%z
        ##     sph.image(s.dm,width=boxsizestring,cmap="Greys", denoise=True ,av_z=False,
        ##               log=True, approximate_fast=False,title=titleStr,
        ##               filename=fileOut); #vmin=0.006, vmax=1.0,                                                                                                    
        ## except:
        ##     print ("EXCEPTION while creating DM image for file: %s"%file)
        ##     gc.collect()
        ##     pass

        # Output Z, PGF, etc. sph images...
        outputImages(z,i,snBlast)

        # Compute stellar masses and save for the mathematica histogram plot
        outputHistogramData(z,i,snBlast)

    del s
    gc.collect()
    gc.collect()

# Sort the keys... Since we'll process the files in Z order
# from low to high, we need to reverse the entries
keys2=key[np.argsort(key[:, 0])]
print ("Keys unsorted: \n",key)
print ("Keys sorted: \n"  ,keys2)
# ##########################################################
np.savetxt("zKeysForSPfiles.txt",keys2,header="z,boxsize_kpc",comments='')
# ##########################################################
# Write out the array to a file
print ("Mass totals\n",data)
np.savetxt("new-data3Mpc.txt",data,comments='',header=headerStr) # comments='' gets rid of "#" in output
del data
print ("Done")
