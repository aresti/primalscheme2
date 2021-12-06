"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains Kmer, Primer and PrimerPair classes.

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

import copy

from collections import namedtuple
from enum import Enum
from typing import Iterable, Sequence

from primalscheme import dna
from primalscheme.config import Config
from primalscheme.interactions import interaction_check


class Kmer:
    """
    0-indexed, always forward, closed end.
    """

    __slots__ = "seq", "start"

    def __init__(self, seq: str, start: int):
        self.seq = seq
        self.start = start

    @property
    def gc(self) -> float:
        return dna.gc_percent(self.seq)

    @property
    def max_homo(self) -> int:
        """Return max homopolymer length for the kmer sequence."""
        return dna.max_homopolymer_length(self.seq)

    def calc_tm(self, cfg: dna.ThermoConfig) -> float:
        """Return Tm for the kmer sequence."""
        return dna.tm(self.seq, cfg)

    def interacts_with(self, kmers: Sequence["Kmer"], cfg: dna.ThermoConfig) -> bool:
        """Return true if the kmer interacts with any of a sequence of kmers."""
        return any(check_kmer_interaction(self, kmer, cfg) for kmer in kmers)

    def passes_thermo_checks(self, cfg: Config) -> bool:
        """Are all kmer thermo values below threshold?"""
        return (
            (cfg.primer_gc_min <= self.gc <= cfg.primer_gc_max)
            and (cfg.primer_tm_min <= self.calc_tm(cfg) <= cfg.primer_tm_max)
            and (self.max_homo <= cfg.primer_homopolymer_max)
        )

    def calc_hairpin_tm(self, cfg: dna.ThermoConfig) -> float:
        """Return the primer3 hairpin thermo object for the kmer sequence."""
        return dna.hairpin(self.seq, cfg).tm

    def passes_hairpin_check(self, cfg: Config) -> bool:
        """Is primer hairpin tm below max?"""
        return self.calc_hairpin_tm(cfg) <= cfg.primer_hairpin_th_max

    def as_reverse(self) -> "Kmer":
        """Return a kmer at the same position with reverse complement seq"""
        kmer = copy.deepcopy(self)
        kmer.seq = self.reverse_complement
        return kmer

    @property
    def reverse_complement(self) -> str:
        """Return the reverse complement fo the seqence."""
        return dna.reverse_complement(self.seq)

    def __len__(self) -> int:
        """Length of kmer sequence."""
        return len(self.seq)

    def __str__(self) -> str:
        return f"Kmer: {self.seq}, start: {self.start}"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Kmer):
            return NotImplemented
        return (self.start, self.seq) == (other.start, other.seq)

    def __hash__(self) -> int:
        return hash((self.start, self.seq))

    @property
    def end(self) -> int:
        """End position (open)."""
        return self.start + len(self)


class PrimerDirection(Enum):
    """Primer direction."""

    FORWARD = "+"
    REVERSE = "-"


class Primer(Kmer):
    """A primer"""

    __slots__ = "direction", "candidate_num"

    def __init__(self, seq: str, start: int, direction: PrimerDirection):
        super().__init__(seq, start)
        self.direction = direction

    @classmethod
    def from_kmer(cls, kmer: Kmer, direction: PrimerDirection) -> "Primer":
        return cls(kmer.seq, kmer.start, direction)

    @classmethod
    def from_bed_row(cls, row: list[str]) -> "Primer":
        seq = row[6]
        direction = (
            PrimerDirection.FORWARD if row[5] == "+" else PrimerDirection.REVERSE
        )
        start = int(row[1])
        return cls(seq, start, direction)

    def __str__(self) -> str:
        return f"{self.direction}, {self.seq}, {self.start}"

    @property
    def end(self) -> int:
        """Primer end position."""
        return self.start + len(self)

    def as_kmer(self) -> Kmer:
        return Kmer(start=self.start, seq=self.seq)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Primer):
            return NotImplemented
        return (self.start, self.seq, self.direction) == (
            other.start,
            other.seq,
            other.direction,
        )

    def __hash__(self) -> int:
        return hash((self.start, self.seq, self.direction))


Insert = namedtuple("Insert", "start end")


class PrimerPair:
    """Primer pair"""

    __slots__ = "forward", "reverse"

    def __init__(self, forward: Primer, reverse: Primer):
        if reverse.start <= forward.end:
            raise ValueError("Reverse primer should be located after forward primer.")

        self.forward = forward
        self.reverse = reverse

    @property
    def start(self) -> int:
        return self.forward.start

    @property
    def end(self) -> int:
        return self.reverse.end

    @property
    def amplicon_size(self) -> int:
        return self.reverse.end - self.forward.start

    @property
    def insert(self) -> Insert:
        return Insert(self.forward.end, self.reverse.start)

    def __len__(self) -> int:
        return self.amplicon_size

    def __str__(self) -> str:
        return f"FWD: {self.forward}, REV: {self.reverse}, SIZE: {self.amplicon_size}"


def digest_seq(seq: str, kmer_size: int) -> Sequence[Kmer]:
    """Digest a sequence into kmers of a given size."""
    return [Kmer(seq[i : i + kmer_size], i) for i in range((len(seq) - kmer_size) + 1)]


def filter_unambiguous_kmers(kmers: Iterable[Kmer]) -> Iterable[Kmer]:
    """Filter: kmers with only unambiguous bases"""
    return [
        kmer
        for kmer in kmers
        if all(base.upper() in dna.UNAMBIGUOUS_DNA for base in kmer.seq)
    ]


def check_kmer_interaction(a: "Kmer", b: "Kmer", cfg: dna.ThermoConfig) -> bool:
    return interaction_check(a.seq, b.seq, cfg) or interaction_check(
        b.reverse_complement, a.reverse_complement, cfg
    )
