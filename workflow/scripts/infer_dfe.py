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
    sfs_sel=spectra['selected']
)

inf.run()

inf.plot_discretized(show=False)
plt.savefig('dfe_discretized.png')

inf.plot_sfs_comparison(show=False)
plt.xlabel('derived allele frequency')
plt.ylabel('number of variants')
plt.tight_layout()
plt.savefig('sfs_comparison.png')

pass
