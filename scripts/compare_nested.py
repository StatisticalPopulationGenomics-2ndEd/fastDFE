"""
Compare nested models.
"""
import fastdfe as fd
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [5, 5]

spectra = fd.Spectra.from_file('sfs.csv')

inf = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    do_bootstrap=True
)

inf_sub = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    do_bootstrap=True,
    fixed_params={'all': {'eps': 0}}
)

p1 = inf_sub.compare_nested(inf)

p2 = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    do_bootstrap=True,
    fixed_params={'all': {'p_b': 0, 'S_b': 1}}
).compare_nested(inf)

# increase margin of plot by default 0.05

inf.plot_nested_models(file='nested_models.png', show=True)

pass
