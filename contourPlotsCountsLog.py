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
def genDensityPlot(x, y, filename, xaxislabel):
    """

    :rtype : none
    """
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111)
    # Note axis order: y then x
    H, xedges, yedges = np.histogram2d(y, x,
                                       range=[[-4., 0.5], [-8, 1.5]], bins=(200, 200))

    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]

    # Use the bins to fine the extent of our plot
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

    plt.subplots_adjust(bottom=0.15, left=0.15)
    levels = (4, 3, 2)

    H = np.log10(H)
    masked_array = np.ma.array(H, mask=np.isnan(H))
    cmap = copy.copy(mpl.cm.jet)
    cmap.set_bad('w', 1.)  # w is color, for values of 1.0

    cset = ax.contour(masked_array, levels, origin='lower', colors=['white', 'black', 'grey'],
                      linewidths=(2.2, 1.9, 1.6), extent=extent, aspect=2.00)

    cax = (ax.imshow(H, extent=extent, cmap=cmap, vmin=0, vmax=4,
                     interpolation='nearest', origin='lower', aspect=2.00))

    cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3, 4])
    cbar.ax.set_yticklabels(['0', '1', '2', '3','4'], size=24)
    cbar.set_label('$Log\, n_{sp}$', size=30)

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
    "05.09",
]
prefix="./"
#prefix="20Sep-BIG/"
for indx, z in enumerate(files):
    spZ = np.loadtxt(prefix+"spZ_" + z + ".txt", skiprows=1)
    spPZ = np.loadtxt(prefix+"spPZ_" + z + ".txt", skiprows=1)
    spPF = np.loadtxt(prefix+"spPPF_" + z + ".txt", skiprows=1)

    print "Generating phase diagram for z=%s"%z
    genDensityPlot(np.log10(spZ), np.log10(spPZ / spZ),
                   "Zplt", "$log\, Z_{\odot}$")
    genDensityPlot(np.log10((spZ) / (1.0-spPF)), np.log10(spPZ / spZ),
                   "Zdiv1-PGFplt", "$log\, Z_{\odot}/f_{pol}$")
