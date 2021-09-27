import math
import pathlib

from typing import Any

UNAMBIGUOUS_DNA = "ACGT"
AMBIGUOUS_DNA = UNAMBIGUOUS_DNA + "RYWSMKHBVDN-"


class Config:
    """
    PrimalDeep configuration.
    Class properties are defaults, can be overriden
    on instantiation (and will shadow class defaults)
    """

    output = pathlib.Path("./output")
    prefix = "scheme"
    force = False

    amplicon_size_min = 380
    amplicon_size_max = 420
    coverage = 2
    packing = 0.5
    high_gc = False

    primer_size_default_min = 19
    primer_size_default_max = 34
    primer_size_default_opt = 22
    primer_size_hgc_min = 17
    primer_size_hgc_max = 30
    primer_size_hgc_opt = 20

    primer_gc_default_min = 30
    primer_gc_default_max = 55
    primer_gc_hgc_min = 40
    primer_gc_hgc_max = 65
    primer_tm_min = 59.5
    primer_tm_max = 62.5
    primer_hairpin_th_max = 47.0
    primer_homopolymer_max = 5

    mv_conc = 100.0
    dv_conc = 2.0
    dntp_conc = 0.8
    dna_conc = 15.0

    def __init__(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

        if self.high_gc:
            self.primer_size_min = self.primer_size_hgc_min
            self.primer_size_max = self.primer_size_hgc_max
            self.primer_size_opt = self.primer_size_hgc_opt
            self.primer_gc_min = self.primer_gc_hgc_min
            self.primer_gc_max = self.primer_gc_hgc_max
        else:
            self.primer_size_min = self.primer_size_default_min
            self.primer_size_max = self.primer_size_default_max
            self.primer_size_opt = self.primer_size_default_opt
            self.primer_gc_min = self.primer_gc_default_min
            self.primer_gc_max = self.primer_gc_default_max

        # Derived
        self.amplicon_size_var = math.ceil(
            (self.amplicon_size_max - self.amplicon_size_min) / 2
        )
        self.insert_size_max = self.amplicon_size_max - (2 * self.primer_size_min)

    def __str__(self) -> str:
        return str(self.__dict__)
