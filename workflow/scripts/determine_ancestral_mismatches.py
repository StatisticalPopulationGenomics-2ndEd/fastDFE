"""
Determine ancestral mismatches between fastDFE's annotation and the pre-defined ancestral alleles.
"""
import ast

import cyvcf2
import fastdfe as fd
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

# get human samples
samples = pd.read_csv("samples.csv")
humans: pd.DataFrame = samples[samples.species == 'homo_sapiens']

f = fd.Filterer(
    vcf="hgdp.anc.vcf.gz",
    filtrations=[
        fd.ExistingOutgroupFiltration(['ref_pan_troglodytes', 'ref_gorilla_gorilla_gorilla', 'ref_pongo_abelii']),
        fd.SNPFiltration(include_samples=humans)
    ],
    output="hgdp.anc.filtered.vcf.gz"
)

f.filter()

reader = cyvcf2.Reader("hgdp.anc.filtered.vcf.gz")

stats = dict(
    mismatches=0,
    matches=0,
    skipped_native=0,
    skipped_comp=0,
    skipped_both=0,
    probs=[],
    probs_mismatches=[],
    info_mismatches=[]
)
for record in tqdm(reader):
    native = record.INFO.get('AA')
    comp = record.INFO.get('AA_ensembl')
    prob = record.INFO.get('AA_prob')
    stats['probs'] += [prob]

    if comp == native and comp:
        stats['matches'] += 1
    elif native == '.' and comp is None:
        stats['skipped_both'] += 1
    elif native == '.':
        stats['skipped_native'] += 1
    elif comp is None:
        stats['skipped_comp'] += 1
    else:
        stats['mismatches'] += 1
        stats['probs_mismatches'] += [prob]
        info = ast.literal_eval(record.INFO.get('AA_info'))
        stats['info_mismatches'] += [info]

plt.hist([p for p in stats['probs'] if p is not None], bins=30, density=True)
plt.xlabel("P(major_ancestral)")
plt.ylabel("Frequency")
plt.title("Major ancestral allele probability distribution")
plt.show()

plt.hist([p for p in stats['probs_mismatches'] if p is not None], bins=30, density=True)
plt.xlabel("P(major_ancestral)")
plt.ylabel("Frequency")
plt.title("Major ancestral allele probability distribution for mismatches")
plt.show()

# calculate mismatch ratio
ratio = stats['mismatches'] / (stats['mismatches'] + stats['matches'])

pass
