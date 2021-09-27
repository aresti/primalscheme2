from collections import namedtuple
from enum import Enum

from .config import Config
from .utils import calc_gc, calc_hairpin, calc_max_homo, calc_tm


"""
Generic feature within the coordinate system
0-based, always forward, closed
"""
Feature = namedtuple("Feature", "start end")


"""
A Kmer
"""
Kmer = namedtuple("Kmer", "seq start")


class PrimerDirection(Enum):
    """Primer direction."""

    FORWARD = "+"
    REVERSE = "-"


class Primer:
    """A primer"""

    def __init__(
        self, seq: str, start: int, direction: PrimerDirection, config: Config
    ) -> None:
        self.cfg = config
        self.seq = seq
        self.start = start
        self.direction = direction

    def __len__(self) -> int:
        """Primer size (length)"""
        return len(self.seq)

    def __str__(self) -> str:
        return f"{self.direction}, {self.seq}, {self.start}"

    @property
    def seq(self) -> str:
        """Primer sequence."""
        return self.__seq

    @seq.setter
    def seq(self, seq: str) -> None:
        """Set seq and calculate derived characteristics."""
        self.__seq = seq
        self.gc = calc_gc(seq)
        self.tm = calc_tm(seq, self.cfg)
        self.hairpin = calc_hairpin(seq, self.cfg).tm
        self.max_homo = calc_max_homo(seq)

    @property
    def end(self) -> int:
        """Primer (inclusive) end position."""
        if self.direction == PrimerDirection.FORWARD:
            return self.start + len(self) - 1
        else:
            return self.start - len(self) + 1


class PrimerPair:
    """Primer pair"""

    def __init__(self, forward: Primer, reverse: Primer, region_center: float):
        self.forward = forward
        self.reverse = reverse
        self.region_center = region_center

    @property
    def start(self) -> int:
        return self.forward.start

    @property
    def end(self) -> int:
        return self.reverse.start

    @property
    def center(self) -> float:
        return self.forward.start + 0.5 * len(self)

    @property
    def position_penalty(self) -> float:
        return abs(self.region_center - self.center)

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __str__(self) -> str:
        return f"{self.forward}, {self.reverse}, {self.position_penalty}"
