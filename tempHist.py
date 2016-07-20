# coding: utf-8
import numpy as np
get_ipython().magic(u'ls totSM*16*.txt')
prefix = "./"
z = "16.0"; i = "1540"
totSPM = np.loadtxt(prefix + "totSM_z" + z + "-" + i +".txt", skiprows=1)
ppfSPM = np.loadtxt(prefix + "pristSM_z" + z + "-" + i +".txt", skiprows=1)
ppzSPM = np.loadtxt(prefix + "primordSM_z" + z + "-" + i +".txt", skiprows=1)
totSPM
ppfSPM
ppzSPM
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import numpy as np
histRect = [0.1, 0.1, .85, 0.85]
axHist = plt.axes(histRect)
N, bins, patches = axHist.hist((totSPM,ppfSPM,ppzSPM), bins=10, log=True)
N
axHist.set_xlim((-10,9))
bins
plt.show()
