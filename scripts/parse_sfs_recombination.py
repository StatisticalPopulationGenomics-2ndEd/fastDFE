"""
Parse SFS with recombination intensity stratification.
"""
import fastdfe as fd
import matplotlib as mpl
import pandas as pd
from matplotlib import pyplot as plt

from scripts.utils import RecombinationIntensityStratification, RecombinationIntensityAnnotation

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
        RecombinationIntensityStratification('recombination_bins.csv')
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
        RecombinationIntensityAnnotation('recombination_map.csv'),
        fd.DegeneracyAnnotation()
    ],
    polarize_probabilistically=True
)

spectra: fd.Spectra = p.parse()
spectra.to_file('sfs.rho.csv')
spectra.plot()

p.annotations[0].plot_histogram(file="rho.hist.png")

spectra.plot(show=False)
plt.legend(ncol=2, fontsize=8)
plt.savefig("sfs.rho.png")

pass
