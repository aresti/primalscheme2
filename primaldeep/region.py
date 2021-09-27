import math

from typing import Iterable, Optional

from .config import Config, UNAMBIGUOUS_DNA
from .exceptions import NoSuitablePrimersError
from .types import Feature, Kmer, Primer, PrimerDirection, PrimerPair
from .utils import primer_thermo_filter, reverse_complement


def digest_seq(seq: str, kmer_size: int) -> list[Kmer]:
    """Digest a sequence into kmers."""
    return [Kmer(seq[i : i + kmer_size], i) for i in range((len(seq) - kmer_size) + 1)]


def filter_unambiguous_kmers(kmers: Iterable[Kmer]) -> Iterable[Kmer]:
    """Filter out kmers with ambiguity codes in their seq"""
    return filter(lambda k: all(b.upper() in UNAMBIGUOUS_DNA for b in k.seq), kmers)


class Region:
    """A defned region within a pool of a scheme"""

    def __init__(
        self,
        seq: str,
        start: int,
        config: Config,
    ):
        self.seq = seq
        self.start = start
        self.cfg = config

        self.fwd_primers: list[Primer] = []
        self.rev_primers: list[Primer] = []
        self.primer_pairs: list[PrimerPair] = []
        self.top_pair: Optional[PrimerPair] = None
        self.center = self.start + 0.5 * len(self.seq)
        self.deviation = 0

        self.run()

    def __str__(self) -> str:
        return f"Center: {self.center}, {self.top_pair}"

    @property
    def fwd_flank_coords(self) -> tuple[int, int]:
        start = math.floor(
            len(self.seq) / 2
            + self.deviation
            - self.cfg.amplicon_size_max / 2
            - self.cfg.primer_size_max
        )
        end = start + self.cfg.amplicon_size_var + self.cfg.primer_size_max
        return start, end

    @property
    def fwd_flank_seq(self) -> str:
        start, end = self.fwd_flank_coords
        return self.seq[start:end]

    @property
    def rev_flank_coords(self) -> tuple[int, int]:
        end = math.ceil(
            len(self.seq) / 2
            + self.deviation
            + self.cfg.amplicon_size_max / 2
            + self.cfg.primer_size_max
        )
        start = end - self.cfg.amplicon_size_var - self.cfg.primer_size_max
        return start, end

    @property
    def rev_flank_seq(self) -> str:
        start, end = self.rev_flank_coords
        return self.seq[start:end]

    @property
    def forward(self) -> Optional[Primer]:
        return self.top_pair.forward if self.top_pair else None

    @property
    def reverse(self) -> Optional[Primer]:
        return self.top_pair.reverse if self.top_pair else None

    @property
    def insert(self) -> Optional[Feature]:
        if not (self.forward and self.reverse):
            return None
        start = self.forward.end + 1
        end = self.reverse.end - 1
        return Feature(start, end)

    def run(self) -> None:
        STEP_DISTANCE = 11
        while True:
            try:
                self.find_primers()
                break
            except NoSuitablePrimersError:
                if self.deviation < 0:
                    self.deviation = abs(self.deviation)
                else:
                    self.deviation = -(self.deviation + STEP_DISTANCE)
                print(f"Stepped to deviation {self.deviation}")

    def find_primers(self) -> None:
        print(f"Flanks: {self.fwd_flank_coords}, {self.rev_flank_coords}")
        print(f"{self.fwd_flank_seq}, {self.rev_flank_seq}")
        self.fwd_primers = []
        self.rev_primers = []
        self.generate_fwd_primers()
        self.generate_rev_primers()
        self.remove_unsuitable_primers()
        print(
            f"Generating pairs from {len(self.fwd_primers)} forward, {len(self.rev_primers)} reverse"
        )
        self.generate_primer_pairs()
        self.sort_primer_pairs()

    def generate_fwd_primers(self) -> None:
        """Generate forward primers for the region"""
        all_kmers = set()

        for kmer_size in range(self.cfg.primer_size_min, self.cfg.primer_size_max + 1):
            all_kmers.update(digest_seq(self.fwd_flank_seq, kmer_size))

        for kmer in filter_unambiguous_kmers(all_kmers):
            fwd_start = self.start + self.fwd_flank_coords[0] + kmer.start
            self.fwd_primers.append(
                Primer(kmer.seq, fwd_start, PrimerDirection.FORWARD, self.cfg)
            )

        if not self.fwd_primers:
            raise NoSuitablePrimersError(
                "No unambiguous kmers within the forward primer zone"
            )

    def generate_rev_primers(self) -> None:
        """Generate all possible reverse primers for the region"""
        all_kmers = set()

        for kmer_size in range(self.cfg.primer_size_min, self.cfg.primer_size_max + 1):
            all_kmers.update(
                digest_seq(reverse_complement(self.rev_flank_seq), kmer_size)
            )

        for kmer in filter_unambiguous_kmers(all_kmers):
            rev_start = self.start + self.rev_flank_coords[1] - 1 - kmer.start
            self.rev_primers.append(
                Primer(kmer.seq, rev_start, PrimerDirection.REVERSE, self.cfg)
            )

        if not self.rev_primers:
            raise NoSuitablePrimersError(
                "No unambiguous kmers within the reverse primer zone"
            )

    def remove_unsuitable_primers(self) -> None:
        self.fwd_primers = list(
            filter(lambda p: primer_thermo_filter(p, self.cfg), self.fwd_primers)
        )
        self.rev_primers = list(
            filter(lambda p: primer_thermo_filter(p, self.cfg), self.rev_primers)
        )
        if not (self.fwd_primers and self.rev_primers):
            raise NoSuitablePrimersError("No suitable primers after thermo filter")

    def generate_primer_pairs(self) -> None:
        """Generate all possible primer pairs for the region"""
        for fwd_primer in self.fwd_primers:
            for rev_primer in self.rev_primers:
                self.primer_pairs.append(
                    PrimerPair(fwd_primer, rev_primer, self.center)
                )

    def sort_primer_pairs(self) -> None:
        self.primer_pairs.sort(key=lambda p: p.position_penalty, reverse=True)
        self.top_pair = self.primer_pairs[0]
