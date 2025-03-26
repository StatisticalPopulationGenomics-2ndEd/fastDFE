"""
Parse unfolded site frequency spectra with degeneracy stratification
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

mpl.rcParams['figure.figsize'] = [8, 3]
mpl.rcParams['figure.dpi'] = 500

# get human samples
samples = pd.read_csv("samples.csv")
humans: pd.DataFrame = samples[samples.species == 'homo_sapiens']

import fastdfe as fd

p = fd.Parser(
    n=20,
    vcf="hgdp.anc.deg.vcf.gz",
    include_samples=humans.id,
    polarize_probabilistically=True,
    stratifications=[fd.DegeneracyStratification()]
)

sfs: fd.Spectra = p.parse()
# sfs.plot()

sfs.plot(show=False)
plt.savefig("sfs.unfolded.deg.png")
sfs.to_file("sfs.unfolded.deg.csv")

pass
