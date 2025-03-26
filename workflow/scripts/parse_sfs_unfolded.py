"""
Parse unfolded SFS.
"""
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [8, 3]
mpl.rcParams['figure.dpi'] = 500

# get human samples
samples = pd.read_csv("samples.csv")
humans: pd.DataFrame = samples[samples.species == 'homo_sapiens']

import fastdfe as fd

p = fd.Parser(
    n=20,
    vcf="hgdp.anc.vcf.gz",
    include_samples=humans.id,
    polarize_probabilistically=True
)

sfs: fd.Spectra = p.parse()
# sfs.plot()

sfs.plot(show=False)
plt.savefig("sfs.unfolded.png")
sfs.to_file('sfs.unfolded.csv')

pass
