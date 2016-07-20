#
# 02 Oct 2015
# Rick Sarmento
#
# Purpose:
#  Reads star particle data and creates contour plots
#  Place histograms of x and y axis along side them
#
# Method:
#
# Revision history
#

# ##########################################################
# Generate colors for histogram bars based on height
# ##########################################################
def colorHistOnHeight(N, patches):
    # we need to normalize the data to 0..1 for the full
    # range of the colormap
    fracs = np.log10(N.astype(float))/8.0 # normalize colors to the top of our scale
    norm = mpl.colors.Normalize(fracs.min(), fracs.max())
    
    for thisfrac, thispatch in zip(fracs, patches):
        color = mpl.cm.jet(thisfrac)
        thispatch.set_facecolor(color)

    return

# ##########################################################
# Generate a combo contour/density plot
# ##########################################################
def genDensityPlot(x, y, mass, pf, filename, xaxislabel):
    """

    :rtype : none
    """
    nullfmt = NullFormatter()

    fig = plt.figure(figsize=(20, 20))
    # axHistx  = plt.subplot2grid((3,3),(0,0), colspan =2)  # Top row, x histogram
    # ax2dhist = plt.subplot2grid((3,3),(1,0), colspan=2, rowspan=2) # 2d histogram
    # axHisty  = plt.subplot2grid((3,3),(1,2), colspan=1, rowspan=2) # Right col, y histogram
    ax2dhist = plt.axes(rect_2dhist)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # Fix any "log10(0)" points...
    x[x == np.inf] = 0.0
    y[y == np.inf] = 0.0
    y[y > 1.0] = 1.0
    # Note axis order: y then x
    H, xedges, yedges = np.histogram2d(y, x, weights=mass * (1.0 - pf),  # Weight counts by polluted mass of sp
                                       range=[[-4., 0.5], [minX, maxX]], bins=(xbins, ybins))

    # print "Edges ", yedges[0], yedges[-1], xedges[0], xedges[-1]
    # Use the bins to fine the extent of our plot
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

    # plt.subplots_adjust(bottom=0.15, left=0.15)
    levels = (5, 4, 3)

    H = np.log10(H)
    masked_array = np.ma.array(H, mask=np.isnan(H))  # mask out all nan, i.e. log10(0.0)
    cmap = copy.copy(mpl.cm.jet)
    cmap.set_bad('w', 1.)  # w is color, for values of 1.0

    # Create the contour plot
    # cset = ax2dhist.contour(masked_array, levels, origin='lower', colors=['red', 'black', 'grey'],
    #                         linewidths=(2.2, 1.9, 1.6), extent=extent, aspect='auto')

    cax = (ax2dhist.imshow(masked_array, extent=extent, cmap=cmap, vmin=0, vmax=8,
                           interpolation='nearest', origin='lower', aspect='auto'))

    cbar = fig.colorbar(cax, ticks=[0, 2, 4, 6, 8])
    cbar.ax.set_yticklabels(['0', '2', '4', '6', '8'], size=24)
    cbar.set_label('$log\, M_{sp, pol,\odot}$', size=30)

    # plt.clabel(cset, inline=1, fontsize=26, fmt='%1.0i')

    # for c in cset.collections:
    #     c.set_linestyle('solid')

    ax2dhist.tick_params(axis='x', labelsize=22)
    ax2dhist.tick_params(axis='y', labelsize=22)
    ax2dhist.set_xlabel(xaxislabel, size=30)
    ax2dhist.set_ylabel('$log\, Z_{pri}/Z$', size=30)
    ax2dhist.set_xlim(ax2dhist.get_xlim())
    ax2dhist.set_ylim(ax2dhist.get_ylim())
    ax2dhist.grid(color='0.75', linestyle=':', linewidth=2)
    
    # Generate the xy axes histograms
    ylims = ax2dhist.get_ylim()
    xlims = ax2dhist.get_xlim()
    ybinList = np.arange(ylims[0],ylims[1],(ylims[1]-ylims[0]+1)/(ybins))
    xbinList = np.arange(xlims[0],xlims[1],(xlims[1]-xlims[0]+1)/(xbins))
    N, bins, patches = axHistx.hist(x, bins=xbinList, log=True, weights=mass * (1.0 - pf))
    colorHistOnHeight(N, patches)
    N, bins, patches = axHisty.hist(y, bins=ybinList, log=True, weights=mass * (1.0 - pf),
                                    orientation='horizontal')
    colorHistOnHeight(N, patches)

    # Setup format of the histograms
    axHistx.set_xlim(ax2dhist.get_xlim())  # Match the x range on the horiz hist
    axHistx.set_ylim([10., 10.e8])  # Constant range for all histograms
    axHistx.tick_params(labelsize=22)
    axHistx.yaxis.set_ticks([1e2,1e4,1e6,1e8])
    axHistx.grid(color='0.75', linestyle=':', linewidth=2)


    axHisty.set_xlim([10., 10.e8])  # We're rotated, so x axis is the value
    axHisty.set_ylim(ax2dhist.get_ylim())  # Match the y range on the vert hist
    axHisty.tick_params(labelsize=22)
    axHisty.xaxis.set_ticks([1e2,1e4,1e6,1e8])
    axHisty.grid(color='0.75', linestyle=':', linewidth=2)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    axHistx.set_title('z=' + z, size=40)

    plt.savefig(filename + "-z_" + z + ".png", dpi=fig.dpi)
    #    plt.show()
    plt.close(fig) # Release memory assoc'd with the plot
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
#import matplotlib.colors as colors # For the colored 1d histogram routine
from matplotlib.ticker import NullFormatter
import numpy as np
import copy as copy

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
# Plot parameters - global
left, width = 0.1, 0.63
bottom, height = 0.1, 0.63
bottom_h = left_h = left + width + 0.01

xbins = ybins = 100

rect_2dhist = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.15]
rect_histy = [left_h, bottom, 0.2, height]

prefix = "./"
for indx, z in enumerate(files):
    spZ = np.loadtxt(prefix + "spZ_" + z + ".txt", skiprows=1)
    spPZ = np.loadtxt(prefix + "spPZ_" + z + ".txt", skiprows=1)
    spPF = np.loadtxt(prefix + "spPPF_" + z + ".txt", skiprows=1)
    spMass = np.loadtxt(prefix + "spMass_" + z + ".txt", skiprows=1)

    print ("Generating phase diagram for z=%s" % z)
    minX = -8.0
    maxX = 0.0
    genDensityPlot(np.log10(spZ), np.log10(spPZ / spZ), spMass, spPF,
                   "Z_pMassZ-MassHistLog", "$log\, Z_{\odot}$")
    minX = -5.0
    genDensityPlot(np.log10((spZ) / (1.0 - spPF)), np.log10(spPZ / spZ), spMass, spPF,
                   "ZmassDiv1-PGFpltMassHist", "$log\, Z_{\odot}/f_{pol}$")
