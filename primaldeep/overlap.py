"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2021 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the OverlapPriorityScheme class.

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

from operator import attrgetter
from typing import Callable, Optional, Sequence

from Bio import SeqRecord  # type: ignore

from primaldeep.exceptions import NoSuitablePrimers
from primaldeep.primer import Kmer, Primer, PrimerDirection, PrimerPair
from primaldeep.scheme import Scheme, ProgressBar
from primaldeep.config import Config


class OverlapPriorityScheme(Scheme):
    def __init__(
        self, ref: SeqRecord, kmers: list[Kmer], cfg: Config, pbar: ProgressBar
    ):
        super().__init__(ref, kmers, cfg, pbar)

        self.pools: tuple[list[PrimerPair], list[PrimerPair]] = ([], [])
        self._pool_num = 0

        self._run()

    @property
    def _pool_index(self) -> int:
        """
        Keeps track of the 'current' pool index while the scheme is running.
        i.e., the index of the pool to which the next PrimerPair will be appened.
        """
        return 1 if len(self.pools[0]) > len(self.pools[1]) else 0

    @property
    def _this_pool(self) -> list[PrimerPair]:
        """Return 'this' pool (to which the next PrimerPair will be appended)."""
        return self.pools[self._pool_index]

    @property
    def _this_pool_primers(self) -> list[Primer]:
        """
        Return a flat list of the primers in this pool.
        """
        return [p for pp in self._this_pool for p in (pp.forward, pp.reverse)]

    @property
    def _other_pool(self) -> list[PrimerPair]:
        """Return the pool that is *not* the current pool."""
        return self.pools[(self._pool_index + 1) % 2]

    @property
    def _prev_pair(self) -> Optional[PrimerPair]:
        """Return the previous pair (the last pair in the 'other' pool)."""
        return self._other_pool[-1] if len(self._other_pool) else None

    @property
    def _prev_pair_same_pool(self) -> Optional[PrimerPair]:
        """
        Return the previous pair in the same pool (last but one pair in the scheme).
        """
        return self._this_pool[-1] if len(self._this_pool) else None

    @property
    def progress(self) -> int:
        """
        Return the end position of prev_pair, representing progress.
        """
        return self._prev_pair.end if self._prev_pair else 0

    def _kmers_right_of_pair(self, pair: PrimerPair) -> Sequence[Kmer]:
        """Return all kmers to the right of a primer pair."""
        return [k for k in self.kmers if k.start > pair.reverse.start]

    def _kmers_maintaining_overlap(
        self, kmers: Sequence[Kmer], pair: PrimerPair
    ) -> Sequence[Kmer]:
        """
        Return all kmers that would maintain min_overlap, if used as forward primer.
        """

        def maintains_overlap(kmer: Kmer) -> bool:
            return (
                kmer.end + self.cfg.min_overlap < pair.reverse.end
                and kmer.start > pair.forward.start
            )

        return [k for k in kmers if maintains_overlap(k)]

    def _overlapping_forward_candidates(self) -> list[Kmer]:
        """
        Return all overlapping, forward candidates, relative to prev_pair,
        sorted by highest end position.
        """
        if self._prev_pair and self._prev_pair_same_pool:
            # Not first in either pool
            non_crashing_kmers = self._kmers_right_of_pair(self._prev_pair_same_pool)
            overlapping_kmers = self._kmers_maintaining_overlap(
                non_crashing_kmers, self._prev_pair
            )
        elif not self._prev_pair:
            # First primer in scheme
            return []
        else:
            # First primer in pool
            overlapping_kmers = self._kmers_maintaining_overlap(
                self.kmers, self._prev_pair
            )
        return sorted(overlapping_kmers, key=attrgetter("end"), reverse=True)

    def _non_overlapping_forward_candidates(self) -> list[Kmer]:
        """Return all non-overlapping, forward candidates, relative to prev_pair."""
        if self._prev_pair:
            return [k for k in self.kmers if k.start >= self._prev_pair.forward.end]
        return self.kmers

    def _forward_candidates(self) -> Sequence[Kmer]:
        """Return all forward candidates (overlapping first, then gapped)"""
        return (
            self._overlapping_forward_candidates()
            + self._non_overlapping_forward_candidates()
        )

    def interaction_checker_factory(self) -> Callable[[Kmer], bool]:
        """
        Return a function that performs the required interaction checks for a candidate
        Kmer against all previously selected primers in the current pool.
        """

        def inner_func(kmer: Kmer) -> bool:
            return not (
                kmer.interacts_with([kmer], self.cfg)
                or kmer.interacts_with(self._this_pool_primers, self.cfg)
            )

        return inner_func

    def _find_next_pair(self) -> PrimerPair:
        """Find the next PrimerPair for the scheme."""
        interaction_checker = self.interaction_checker_factory()

        for fwd in self._forward_candidates():
            if interaction_checker(fwd):
                candidate_pairs = self.reverse_candidate_pairs(
                    Primer(fwd.seq, fwd.start, PrimerDirection.FORWARD)
                )
                for pair in candidate_pairs:
                    if interaction_checker(pair.reverse):
                        return pair

        raise NoSuitablePrimers

    def _run(self) -> None:
        """Create a an overlap-priority scheme."""

        last_progress = 0

        while True:
            try:
                self._this_pool.append(self._find_next_pair())
                self.pbar.update(self.progress - last_progress)
                last_progress = self.progress
            except NoSuitablePrimers:
                break
