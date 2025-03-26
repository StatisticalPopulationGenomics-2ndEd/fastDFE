"""
Annotate the ancestral allele of human samples in the HGDP dataset.
"""
import fastdfe as fd
import pandas as pd

# get human samples
samples = pd.read_csv("samples.csv")
humans: pd.DataFrame = samples[samples.species == 'homo_sapiens']

anc = fd.MaximumLikelihoodAncestralAnnotation(
    model=fd.K2SubstitutionModel(
        fix_transition_transversion_ratio=True
    ),
    outgroups=[
        'ref_pan_troglodytes',
        'ref_gorilla_gorilla_gorilla',
        'ref_pongo_abelii'
    ],
    ingroups=humans,
    max_sites=7500,
    n_target_sites=3e6,
    prior=fd.AdaptivePolarizationPrior()
)

ann = fd.Annotator(
    annotations=[anc],
    vcf="hgdp.vcf.gz",
    fasta="hg38.fasta.gz",
    output="hgdp.anc.vcf.gz",
)

ann.annotate()

# anc.prior.plot()

anc.to_file("anc.json")

r = anc.get_observed_transition_transversion_ratio()
div = anc.get_outgroup_divergence()

pass
