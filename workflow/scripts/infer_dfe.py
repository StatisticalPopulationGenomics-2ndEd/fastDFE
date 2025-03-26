"""
Infer the DFE.
"""
import fastdfe as fd
import matplotlib as mpl
from matplotlib import pyplot as plt

mpl.rcParams['figure.figsize'] = [8, 3]
mpl.rcParams['figure.dpi'] = 500

spectra = fd.Spectra.from_file('sfs.csv')

inf = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    do_bootstrap=True,
    n_bootstraps=100
)

inf.run()

inf.plot_discretized(show=False)
plt.savefig('dfe_discretized.png')

inf.plot_sfs_comparison(show=False)
plt.savefig('sfs_comparison.png')

pass
