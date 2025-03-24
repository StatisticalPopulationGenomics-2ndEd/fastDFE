"""
Parse the SFS using random stratification.
This is done to get an idea on the sampling variance of the SFS.
"""
import fastdfe as fd
import matplotlib as mpl
import pandas as pd
from matplotlib import pyplot as plt

mpl.rcParams['figure.figsize'] = [8, 3]

# get human samples
samples = pd.read_csv("samples.csv")
humans: pd.DataFrame = samples[samples.species == 'homo_sapiens']

p = fd.Parser(
    n=20,
    vcf="hgdp.anc.deg.rho.vcf.gz",
    fasta="hg38.fasta.gz",
    gff="hg38.gff3.gz",
    include_samples=humans.id,
    stratifications=[
        fd.DegeneracyStratification(),
        fd.RandomStratification(n_bins=9)
    ],
    filtrations=[
        fd.SNPFiltration(),
        fd.ExistingOutgroupFiltration([
            'ref_pan_troglodytes',
            'ref_gorilla_gorilla_gorilla',
            'ref_pongo_abelii'
        ])
    ],
    target_site_counter=fd.TargetSiteCounter(
        n_target_sites=2e7,
        n_samples=1e8
    ),
    annotations=[
        fd.DegeneracyAnnotation()
    ],
    polarize_probabilistically=True
)

spectra: fd.Spectra = p.parse()
spectra.to_file('sfs.random.csv')
spectra.plot()

spectra.plot(show=False)
plt.legend(ncol=2, fontsize=8)
plt.savefig("sfs.random.png")

pass
