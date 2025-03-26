"""
Parse SFS for chromosome 1 and 2.
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
    fasta="hg38.fasta.gz",
    gff="hg38.gff3.gz",
    include_samples=humans.id,
    stratifications=[
        fd.DegeneracyStratification(),
        fd.ContigStratification(['chr1', 'chr2'])
    ],
    filtrations=[
        fd.ContigFiltration(['chr1', 'chr2']),
        fd.SNPFiltration(),
        fd.ExistingOutgroupFiltration([
            'ref_pan_troglodytes',
            'ref_gorilla_gorilla_gorilla',
            'ref_pongo_abelii'
        ])
    ],
    target_site_counter=fd.TargetSiteCounter(
        n_target_sites=2e7 / 6,
        n_samples=2e6
    ),
    annotations=[fd.DegeneracyAnnotation()],
    polarize_probabilistically=True
)

sfs: fd.Spectra = p.parse()
# sfs.plot()
sfs.to_file('sfs.chr.csv')

sfs.plot(show=False)
plt.savefig("sfs.chr.png")

pass
