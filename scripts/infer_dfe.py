"""
Infer the DFE.
"""
import fastdfe as fd
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [8, 3]

spectra = fd.Spectra.from_file('sfs.csv')

inf = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    do_bootstrap=True,
    n_bootstraps=100
)

inf.run()

inf.plot_discretized(file='dfe_discretized.png')
inf.plot_sfs_comparison(file='sfs_comparison.png')

pass
