"""
Parse SFS including monomorphic site counts.
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
    vcf="hgdp.anc.deg.vcf.gz",
    fasta="hg38.fasta.gz",
    gff="hg38.gff3.gz",
    include_samples=humans.id,
    stratifications=[fd.DegeneracyStratification()],
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
        n_samples=1e6
    ),
    annotations=[fd.DegeneracyAnnotation()],
    polarize_probabilistically=True
)

sfs: fd.Spectra = p.parse()
# sfs.plot()

sfs.plot(show=False)
plt.savefig("sfs.png")
sfs.to_file('sfs.csv')

pass
