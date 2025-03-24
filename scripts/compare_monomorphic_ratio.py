"""
Compare the monomorphic ratio of selected and neutral regions.
"""
import fastdfe as fd

spectra = fd.Spectra.from_file('sfs.rho.csv')

df = spectra._to_multi_index().data
ratios = df.xs('selected', axis=1, level=0).iloc[0] / df.xs('neutral', axis=1, level=0).iloc[0]

r = ratios.std()

pass
