"""
Parse unstratified SFS.
"""
import pandas as pd
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8, 3]

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
sfs.fold().plot()  # plot folded spectrum

sfs.fold().plot(show=False, file="sfs_plain.png")
sfs.to_file('sfs_plain.csv')

pass
