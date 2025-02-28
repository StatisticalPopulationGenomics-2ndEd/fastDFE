"""
Annotation of recombination intensity.
"""
import fastdfe as fd

from workflow.scripts.utils import RecombinationIntensityAnnotation

ann = fd.Annotator(
    annotations=[RecombinationIntensityAnnotation('recombination_map.csv', 'recombination_bins.csv')],
    vcf="hgdp.anc.deg.vcf.gz",
    gff="hg38.gff3.gz",
    fasta="hg38.fasta.gz",
    output="hgdp.anc.deg.rho.vcf.gz"
)

ann.annotate()

pass
