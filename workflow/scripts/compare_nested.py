"""
Compare nested models.
"""
import fastdfe as fd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['figure.figsize'] = [5, 5]
mpl.rcParams['figure.dpi'] = 500

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
print(f'p-value: {p1}')

p2 = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    do_bootstrap=True,
    fixed_params={'all': {'p_b': 0, 'S_b': 1}}
).compare_nested(inf)
print(f'p-value: {p2}')

# increase margin of plot by default 0.05

inf.plot_nested_models(show=False)
plt.savefig('nested_models.png')

pass
