"""
Infer joint DFE.
"""
import fastdfe as fd
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [8, 3]

spectra = fd.Spectra.from_file('sfs.chr.csv')

inf = fd.JointInference(
    sfs_neut=spectra['neutral.*'].rename(['chr1', 'chr2']),
    sfs_sel=spectra['selected.*'].rename(['chr1', 'chr2']),
    fixed_params={'all': {'eps': 0, 'p_b': 0, 'S_b': 1}},
    shared_params=[fd.SharedParams(types=['chr1', 'chr2'], params=['S_d'])],
    do_bootstrap=True,
    n_bootstraps=100
)

inf.run()
inf.plot_discretized(file='dfe_discretized.joint.png')

p = inf.perform_lrt_shared()

pass
