"""
Parse unfolded site frequency spectra with degeneracy stratification
"""
import pandas as pd
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8, 3]

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
sfs.plot()

sfs.plot(show=False, file="sfs.unfolded.deg.png")
sfs.to_file("sfs.unfolded.deg.csv")

pass