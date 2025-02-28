"""
Showcase the different DFE parametrizations.
"""
import fastdfe as fd
import matplotlib as mpl
import numpy as np

mpl.rcParams['figure.figsize'] = [8, 3]
mpl.rcParams['legend.loc'] = 'upper right'

spectra = fd.Spectra.from_file('sfs.csv')

parametrizations = [
    fd.GammaExpParametrization(),
    fd.DiscreteFractionalParametrization(),
    fd.GammaDiscreteParametrization(),
    fd.DisplacedGammaParametrization()
]

inferences = [fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    do_bootstrap=True,
    model=p
) for p in parametrizations]

[inf.run() for inf in inferences]

fd.Inference.plot_discretized(
    inferences=inferences,
    labels=[type(p).__name__ for p in parametrizations],
    intervals=[-np.inf, -100, -10, -1, 1, np.inf],
    file='dfe_parametrizations.png'
)

pass
