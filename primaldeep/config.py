"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2021 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the Config class.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"""

import statistics
import pathlib

from typing import Any, Optional


class Config:
    """
    PrimalDeep configuration.
    Class properties are defaults, can be overriden
    on instantiation (and will shadow class defaults)
    """

    output = pathlib.Path("./output")
    prefix = "scheme"
    force = False
    repair: Optional[pathlib.Path] = None

    amplicon_size_min = 380
    amplicon_size_max = 420
    amplicon_size_target: int
    min_overlap = 10
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
    dimer_max_tm = -10.0
    dimer_min_identity = 0.8

    def __init__(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

        self.amplicon_size_target = int(
            statistics.mean([self.amplicon_size_min, self.amplicon_size_max])
        )

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

    def __str__(self) -> str:
        return str(self.__dict__)
