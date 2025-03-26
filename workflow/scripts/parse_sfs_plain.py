"""
Parse unstratified SFS.
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
    n=20,  # number of SFS bins
    vcf="hgdp.vcf.gz",
    skip_non_polarized=False,  # allow non-polarized sites
    include_samples=humans.id  # only include human samples
)

sfs: fd.Spectra = p.parse()
# sfs.fold().plot()  # plot folded spectrum

sfs.fold().plot(show=False)
plt.savefig("sfs.plain.png")
sfs.to_file('sfs.plain.csv')

pass
