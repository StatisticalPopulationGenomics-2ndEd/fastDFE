"""
Plot an SFS.
"""
import fastdfe as fd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['figure.figsize'] = [8, 3.2]
mpl.rcParams['figure.dpi'] = 500

sfs = fd.Spectra.from_file('sfs.rho.csv')
sfs.plot(show=False)
plt.xlabel('minor allele frequency')
plt.ylabel('number of variants')
plt.legend(ncol=2, fontsize=8)
plt.tight_layout(pad=2)
plt.savefig("sfs.rho.png")

pass
