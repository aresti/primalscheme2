"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2021 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains utilities for dealing with DNA strings, and wrapped primer
functions.

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

import re

from itertools import groupby
from typing_extensions import Protocol

from primer3 import calcTm as p3_calc_tm, calcHairpin as p3_calc_hairpin  # type: ignore

UNAMBIGUOUS_DNA = "ACGT"
AMBIGUOUS_DNA = UNAMBIGUOUS_DNA + "RYWSMKHBVDN-"
AMBIGUOUS_DNA_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
}

CIGAR_REGEX = re.compile(r"(?P<len>\d+)(?P<op>\D+)")


class SeqProtocol(Protocol):
    seq: str

    def __len__(self) -> int:
        ...


class SeqRecordProtocol(Protocol):
    id: str
    seq: SeqProtocol

    def format(self, format: str) -> str:
        ...


class ThermoConfig(Protocol):
    mv_conc: float
    dv_conc: float
    dntp_conc: float
    dna_conc: float
    dimer_max_tm: float
    dimer_min_identity: float


class ThermoResult(Protocol):
    structure_found: bool
    tm: float
    dg: float
    dh: float
    ds: float


def gc_percent(seq: str) -> float:
    """Caclulate percent GC for a sequence."""
    return round(100.0 * (seq.count("G") + seq.count("C")) / len(seq), 1)


def tm(seq: str, config: ThermoConfig) -> float:
    """Calculate Tm for a sequence using primer3."""
    return p3_calc_tm(
        seq,
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dntp_conc=config.dntp_conc,
        dna_conc=config.dna_conc,
    )


def hairpin(seq: str, config: ThermoConfig) -> ThermoResult:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Return primer3 ThermoResult object.
    """
    return p3_calc_hairpin(
        seq,
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dntp_conc=config.dntp_conc,
        dna_conc=config.dna_conc,
    )


def max_homopolymer_length(seq: str) -> int:
    """Calculate max homopolymer length for a sequence."""
    if not seq:
        return 0
    return max(sum(1 for _ in group) for _, group in groupby(seq))


def reverse(seq: str) -> str:
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq: str) -> str:
    """Returns the complement of seq"""
    if not all(base.upper() in UNAMBIGUOUS_DNA for base in seq):
        raise ValueError(f"Sequence bases must be one of {UNAMBIGUOUS_DNA}")
    return "".join(AMBIGUOUS_DNA_COMPLEMENT[base.upper()] for base in seq)


def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of seq"""
    return reverse(complement(seq))
