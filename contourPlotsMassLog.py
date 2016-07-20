#
# 02 Oct 2015
# Rick Sarmento
#
# Purpose:
#  Reads star particle data and creates contour plots
#
# Method:
#
# Revision history
#

# ##########################################################
# Generate a combo contour/density plot
# ##########################################################
def genDensityPlot(x, y, mass, pf, filename, xaxislabel):
    """

    :rtype : none
    """
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111)
    # Note axis order: y then x
    H, xedges, yedges = np.histogram2d(y, x, weights= mass * (1.0 - pf),
                                       range=[[-4., 0.5], [-8, 1.5]], bins=(200, 200))

    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]

    #print "Edges ", yedges[0], yedges[-1], xedges[0], xedges[-1]
    # Use the bins to fine the extent of our plot
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

    plt.subplots_adjust(bottom=0.15, left=0.15)
    levels = (6, 4, 3)

    H = np.log10(H)
    masked_array = np.ma.array(H, mask=np.isnan(H))
    cmap = copy.copy(mpl.cm.jet)
    cmap.set_bad('w', 1.)  # w is color, for values of 1.0

    cset = ax.contour(masked_array, levels, origin='lower', colors=['white', 'black', 'grey'],
                      linewidths=(2.2, 1.9, 1.6), extent=extent, aspect=2.00)

    cax = (ax.imshow(H, extent=extent, cmap=cmap, vmin=0, vmax=12,
                     interpolation='nearest', origin='lower', aspect=2.00))

    cbar = fig.colorbar(cax, ticks=[0, 2, 4, 6, 8])
    cbar.ax.set_yticklabels(['0', '2', '4', '6','8'], size=24)
    cbar.set_label('$Log\, m_{sp, pol}$', size=30)

    plt.clabel(cset, inline=1, fontsize=26, fmt='%1.0i')

    for c in cset.collections:
        c.set_linestyle('solid')

    ax.tick_params(axis='x', labelsize=24)
    ax.tick_params(axis='y', labelsize=24)
    ax.set_title('z=' + z, size=40)
    ax.set_xlabel(xaxislabel, size=30)
    ax.set_ylabel('$log\, Z_{pri}/Z$', size=30)
    plt.savefig(filename+"-z_" + z + ".png", dpi=fig.dpi)
    #    plt.show()

    plt.close()
    return


# ##########################################################
# ##########################################################
##
## Main program
##
# ##########################################################
# ##########################################################
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import copy

files = [
    "18.00",
    "17.00",
    "16.00",
    "15.00",
    "14.00",
    "13.00",
    "12.00",
    "11.00",
    "10.00",
    "09.00",
    "08.50",
    "08.00",
    "07.50",
    "07.00",
    "06.50",
    "06.00",
    "05.50",
    "05.09"
]
prefix="./"
#prefix="20Sep-BIG/"
for indx, z in enumerate(files):
    spZ = np.loadtxt(prefix+"spZ_" + z + ".txt", skiprows=1)
    spPZ = np.loadtxt(prefix+"spPZ_" + z + ".txt", skiprows=1)
    spPF = np.loadtxt(prefix+"spPPF_" + z + ".txt", skiprows=1)
    spMass = np.loadtxt(prefix+"spMass_" + z + ".txt", skiprows=1)

    print "Generating phase diagram for z=%s"%z
    print "largest mass at this z is %.1lf"%max(spMass)
    genDensityPlot(np.log10(spZ), np.log10(spPZ / spZ), spMass, spPF,
                   "ZmassPlt", "$log\, Z_{\odot}$")
    genDensityPlot(np.log10((spZ) / (1.0-spPF)), np.log10(spPZ / spZ), spMass, spPF,
                   "ZmassDiv1-PGFplt", "$log\, Z_{\odot}/f_{pol}$")
