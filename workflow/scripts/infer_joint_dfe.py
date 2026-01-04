"""
Infer joint DFE.
"""
import fastdfe as fd
import matplotlib as mpl
from matplotlib import pyplot as plt

mpl.rcParams['figure.figsize'] = [8, 3]
mpl.rcParams['figure.dpi'] = 500

spectra = fd.Spectra.from_file('sfs.chr.csv')

# deleterious DFE, shared S_d
inf = fd.JointInference(
    sfs_neut=spectra['neutral.*'].rename(['chr1', 'chr2']),
    sfs_sel=spectra['selected.*'].rename(['chr1', 'chr2']),
    fixed_params={'all': {'h': 0.5, 'eps': 0, 'p_b': 0, 'S_b': 1}},
    shared_params=[fd.SharedParams(types=['chr1', 'chr2'], params=['S_d'])]
)

inf.run()
inf.plot_discretized(show=False)
plt.savefig('dfe_discretized.joint.png')

p = inf.perform_lrt_shared()
print(f'p-value (deleterious DFE, shared S_d): {p}')

# full DFE, shared S_d
inf = fd.JointInference(
    sfs_neut=spectra['neutral.*'].rename(['chr1', 'chr2']),
    sfs_sel=spectra['selected.*'].rename(['chr1', 'chr2']),
    fixed_params={'all': {'h': 0.5, 'eps': 0}},
    shared_params=[fd.SharedParams(types=['chr1', 'chr2'], params=['S_d'])]
)

inf.run()
p1 = inf.perform_lrt_shared()
print(f'p-value (full DFE, shared S_d): {p1}')

# deleterious DFE, shared b
inf = fd.JointInference(
    sfs_neut=spectra['neutral.*'].rename(['chr1', 'chr2']),
    sfs_sel=spectra['selected.*'].rename(['chr1', 'chr2']),
    fixed_params={'all': {'h': 0.5, 'eps': 0, 'p_b': 0, 'S_b': 1}},
    shared_params=[fd.SharedParams(types=['chr1', 'chr2'], params=['b'])]
)

inf.run()
p2 = inf.perform_lrt_shared()
print(f'p-value (deleterious DFE, shared b): {p2}')

pass
