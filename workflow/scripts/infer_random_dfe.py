"""
Infer the DFE with covariates stratified randomly.
"""
import fastdfe as fd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams['figure.figsize'] = [8, 4]
mpl.rcParams['figure.dpi'] = 500

spectra = fd.Spectra.from_file('sfs.random.csv')
bins = pd.read_csv('rho_bins.csv').iloc[:, 1].to_numpy()
values = np.log((bins[:-1] + bins[1:]) / 2)
types = spectra.merge_groups(level=1).types
covariate = fd.Covariate('S_d', dict(zip(types, values)))

pass

inf = fd.JointInference(
    sfs_neut=spectra['neutral.*'].merge_groups(level=1),
    sfs_sel=spectra['selected.*'].merge_groups(level=1),
    fixed_params={'all': {'eps': 0, 'p_b': 0, 'S_b': 1}},
    shared_params=[fd.SharedParams(['b', 'S_d'])],
    covariates=[covariate],
    do_bootstrap=True,
    n_bootstraps=100
)

inf.run()
inf.plot_discretized(show=False)
plt.legend(ncol=2, fontsize=7)
plt.savefig('dfe_discretized.random.png')

c0 = inf.bootstraps.mean()['bin0.c0']
print(f'c0: {c0}')

S_d = {k: v.bootstraps.mean()['S_d'] for k, v in inf.marginal_inferences.items() if k != 'all'}
plt.bar(S_d.keys(), np.abs(list(S_d.values())))
plt.ylabel('$S_d$')
plt.xlabel('Recombination intensity bin')
plt.title('Marginal inferences')
plt.tight_layout()
plt.savefig('marginal_inferences.random.png')

S_d = {k: v.bootstraps.mean()['S_d'] for k, v in inf.joint_inferences.items() if k != 'all'}
plt.bar(S_d.keys(), np.abs(list(S_d.values())))
plt.ylabel('$S_d$')
plt.xlabel('Recombination intensity bin')
plt.title('Joint inferences')
plt.tight_layout()
plt.savefig('joint_inferences.random.png')

mpl.rcParams['figure.figsize'] = [8, 3]
p = inf.perform_lrt_covariates()
inf.plot_covariate(show=False)
plt.savefig('lrt_covariates.random.png')

pass
