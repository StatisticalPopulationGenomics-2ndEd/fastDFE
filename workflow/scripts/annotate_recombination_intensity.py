"""
Annotation of recombination intensity.
"""
import os

import fastdfe as fd
import sys

sys.path.append(os.getcwd() + "/..")

from workflow.scripts.utils import RecombinationIntensityAnnotation

ann = fd.Annotator(
    annotations=[RecombinationIntensityAnnotation('rho_map.csv', 'rho_bins.csv')],
    vcf="hgdp.anc.deg.vcf.gz",
    gff="hg38.gff3.gz",
    fasta="hg38.fasta.gz",
    output="hgdp.anc.deg.rho.vcf.gz"
)

ann.annotate()

pass
