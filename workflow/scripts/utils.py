"""
Utility classes.
"""
from typing import Union, List

import fastdfe as fd
import numpy as np
import pandas as pd
from fastdfe.io_handlers import DummyVariant


class RecombinationIntensityAnnotation(fd.Annotation):
    """
    Recombination intensity annotation using a genetic map.
    """

    def __init__(self, map_file: str, bins_file: str = None, n_bins: int = 10):
        """
        Create a new annotation instance.
        """
        super().__init__()

        self.bins_file = bins_file

        self.n_bins = n_bins

        self.map: pd.DataFrame = pd.read_csv(map_file, comment='#', sep='\t')

        self.chr_cache = {chr_: df for chr_, df in self.map.groupby("Chr")}

        self.values: List[float] = []

        self.current: pd.Series = self.map.iloc[0]

    def _setup(self, handler: fd.io_handlers.MultiHandler):
        """
        Provide context to the annotator.
        """
        super()._setup(handler)

        handler._reader.add_info_to_header({
            'ID': 'Recombination_intensity',
            'Number': '.',
            'Type': 'Float',
            'Description': 'Recombination intensity in cM/Mb.'
        })

    def _teardown(self):
        """
        Finalize the annotation. Called after all sites have been annotated.
        """
        super()._teardown()

        if self.bins_file is not None:
            pd.DataFrame(np.quantile(self.values, np.linspace(0, 1, self.n_bins + 1))).to_csv(self.bins_file)

    def annotate_site(self, variant: Union['cyvcf2.Variant', DummyVariant]):
        """
        Annotate site.
        """
        if (self.current is not None and variant.CHROM == self.current.Chr and
                self.current.Begin <= variant.POS < self.current.End):
            variant.INFO['Recombination_intensity'] = float(self.current.cMperMb)
            self.values.append(self.current.cMperMb)
            self.n_annotated += 1
            return

        # get chromosome-specific dataframe
        chr_map = self.chr_cache.get(variant.CHROM)
        if chr_map is not None:
            mask = (chr_map.Begin <= variant.POS) & (variant.POS < chr_map.End)
            if mask.any():
                self.current = chr_map[mask].iloc[0]
                self.annotate_site(variant)


class RecombinationIntensityStratification(fd.Stratification):
    """
    Recombination intensity stratification based annotation.
    """

    def __init__(self, bins_file: str):
        """
        Initialize stratification.
        """
        super().__init__()

        self.bins: pd.Series = pd.read_csv(bins_file).iloc[:, 1]

    def get_types(self) -> List[str]:
        """
        Get all possible types.
        """
        return [f"bin{i}" for i in range(len(self.bins) - 1)]

    def get_type(self, variant: Union['cyvcf2.Variant', DummyVariant]) -> str:
        """
        Get type.
        """
        try:
            rho = variant.INFO['Recombination_intensity']
            return f"bin{sum(rho >= self.bins) - 1}"
        except KeyError:
            raise fd.io_handlers.NoTypeException()
