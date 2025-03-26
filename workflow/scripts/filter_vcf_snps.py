"""
Filter VCF file to keep only SNPs and human samples.
"""
import fastdfe as fd
import pandas as pd

# get human samples
samples = pd.read_csv("samples.csv")
humans: pd.DataFrame = samples[samples.species == 'homo_sapiens']

f = fd.Filterer(
    vcf="hgdp.vcf.gz",
    filtrations=[
        fd.SNPFiltration(include_samples=humans)
    ],
    output="hgdp.snps.vcf.gz"
)

f.filter()

pass
