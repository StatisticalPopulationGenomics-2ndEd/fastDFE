"""
Annotation of degeneracy for the HGDP dataset.
"""
import fastdfe as fd

ann = fd.Annotator(
    annotations=[fd.DegeneracyAnnotation()],
    vcf="hgdp.anc.vcf.gz",
    gff="hg38.gff3.gz",
    fasta="hg38.fasta.gz",
    output="hgdp.anc.deg.vcf.gz"
)

ann.annotate()
