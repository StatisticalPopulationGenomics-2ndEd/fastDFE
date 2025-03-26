"""
Filter VCF file to retain only 0 and 4-fold degenerate sites.
"""
import fastdfe as fd
from cyvcf2 import Reader, Writer
from tqdm import tqdm

ann = fd.Annotator(
    annotations=[fd.DegeneracyAnnotation()],
    vcf="hgdp.raw.vcf.gz",
    gff="hg38.gff3.gz",
    fasta="hg38.fasta.gz",
    output="hgdp.raw.deg.vcf.gz"
)

ann.annotate()

# only retain 0 and 4-fold degenerate sites and remove degeneracy tag
reader = Reader("hgdp.raw.deg.vcf.gz")
writer = Writer("hgdp.vcf.gz", reader)

# Filter and write variants
for variant in tqdm(reader):
    if variant.INFO.get('Degeneracy') in [0, 2, 4]:
        del variant.INFO['Degeneracy']
        del variant.INFO['Degeneracy_Info']
        writer.write_record(variant)

# Close writer and reader
writer.close()
reader.close()

pass
