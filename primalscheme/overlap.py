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
from loguru import logger
from operator import attrgetter
from typing import Callable, Optional, Sequence

from Bio import SeqRecord  # type: ignore

from primalscheme.exceptions import NoSuitablePrimers
from primalscheme.primer import Kmer, Primer, PrimerDirection, PrimerPair
from primalscheme.scheme import Scheme
from primalscheme.config import Config

logger = logger.opt(colors=True)


class OverlapPriorityScheme(Scheme):
    def __init__(
        self, ref: SeqRecord, fwd_kmers: list[Kmer], rev_kmers: list[Kmer], cfg: Config
    ):
        super().__init__(ref, fwd_kmers, rev_kmers, cfg)

        self.pools: tuple[list[PrimerPair], list[PrimerPair]] = ([], [])
        self._pool_num = 0

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

    def _available_fwd_kmers(self) -> Sequence[Kmer]:
        """Return all available forward kmers."""
        if self._prev_pair is None:
            # Empty scheme
            return self.fwd_kmers

        if self._prev_pair_same_pool is None:
            # Empty pool
            start_from = self._prev_pair.forward.start + 1
        else:
            # Not first in either pool
            start_from = self._prev_pair_same_pool.reverse.end
        return [k for k in self.fwd_kmers if k.start >= start_from]

    def _fwd_kmers_maintaining_overlap(
        self,
    ) -> list[Kmer]:
        """
        Return all kmers that would maintain min_overlap, if used as forward primer,
        sorted by highest end position (furthest right).
        """
        if self._prev_pair is None:
            return []

        pair = self._prev_pair

        def maintains_overlap(kmer: Kmer) -> bool:
            return (
                kmer.end - 1 + self.cfg.min_overlap < pair.reverse.start
                and kmer.start > pair.forward.start
            )

        overlapping = [k for k in self._available_fwd_kmers() if maintains_overlap(k)]
        return sorted(overlapping, key=attrgetter("end"), reverse=True)

    def _fwd_kmers_non_overlapping(self) -> list[Kmer]:
        """Return all non-overlapping, forward candidates, relative to prev_pair."""
        if self._prev_pair is None:
            return self.fwd_kmers
        return [k for k in self.fwd_kmers if k.start >= self._prev_pair.forward.end]

    def _fwd_candidates(self) -> list[Kmer]:
        """Return all forward candidates (overlapping first, then non-overlapping)"""
        return self._fwd_kmers_maintaining_overlap() + self._fwd_kmers_non_overlapping()

    def interaction_checker_factory(
        self, existing_pairs: Optional[list[PrimerPair]] = None
    ) -> Callable[[Kmer], bool]:
        """
        Return a function that performs the required interaction checks for a candidate
        Kmer against all previously selected primers in the current pool.
        """

        check_primers = self._this_pool_primers
        if existing_pairs is not None:
            check_primers += [
                p for pp in existing_pairs for p in (pp.forward, pp.reverse)
            ]

        def inner_func(kmer: Kmer) -> bool:
            return not (
                kmer.interacts_with([kmer], self.cfg)
                or kmer.interacts_with(check_primers, self.cfg)
            )

        return inner_func

    def _find_pair(
        self,
        next_pair: Optional[PrimerPair] = None,
        same_pool_pairs: Optional[list[PrimerPair]] = None,
    ) -> PrimerPair:
        """Find the next PrimerPair for the scheme."""
        interaction_checker = self.interaction_checker_factory(
            existing_pairs=same_pool_pairs
        )

        for fwd in self._fwd_candidates():
            if interaction_checker(fwd):
                try:
                    candidate_pairs = self.reverse_candidate_pairs(
                        Primer.from_kmer(fwd, PrimerDirection.FORWARD),
                        next_pair=next_pair,
                    )
                except NoSuitablePrimers:
                    continue
                for pair in candidate_pairs:
                    if interaction_checker(pair.reverse):
                        return pair

        raise NoSuitablePrimers

    def execute(self) -> None:
        """Create a an overlap-priority scheme."""
        logger.info("Designing scheme with <blue>overlap-priority</> strategy")
        while True:
            try:
                self._this_pool.append(self._find_pair())
            except NoSuitablePrimers:
                break

    def _replacement_fwd_kmers(self, pair: PrimerPair) -> list[Kmer]:
        """
        Return replacement forward candidate Kmers for an existing pair.
        """
        reverse = pair.reverse
        window = (
            reverse.end - self.cfg.amplicon_size_max - 1,
            reverse.end - self.cfg.amplicon_size_min,
        )
        return [
            k
            for k in self._fwd_kmers_maintaining_overlap()
            if k.start >= window[0] and k.end <= window[1]
        ]

    def _replacement_fwd_pairs(self, pair: PrimerPair) -> list[PrimerPair]:
        rev = pair.reverse
        pairs = [
            PrimerPair(Primer.from_kmer(kmer, PrimerDirection.FORWARD), rev)
            for kmer in self._replacement_fwd_kmers(pair)
        ]

        # Sort by deviation from target amplicon size
        pairs.sort(key=lambda p: abs(self.cfg.amplicon_size_target - p.amplicon_size))
        return pairs

    def _find_replacement_fwd(
        self, pair: PrimerPair, existing_pairs: list[PrimerPair]
    ) -> PrimerPair:
        interaction_check = self.interaction_checker_factory(
            existing_pairs=existing_pairs
        )
        for candidate in self._replacement_fwd_pairs(pair):
            if interaction_check(candidate.forward):
                return candidate
        raise NoSuitablePrimers

    def _replacement_rev_pairs(
        self, pair: PrimerPair, next_pair: PrimerPair
    ) -> list[PrimerPair]:
        fwd = pair.forward
        return [
            pp
            for pp in self.reverse_candidate_pairs(fwd)
            if pp.reverse.start >= next_pair.forward.end
        ]

    def _find_replacement_rev(
        self,
        pair: PrimerPair,
        next_pair: Optional[PrimerPair] = None,
        existing_pairs: Optional[list[PrimerPair]] = None,
    ) -> PrimerPair:
        interaction_check = self.interaction_checker_factory(
            existing_pairs=existing_pairs
        )
        if next_pair is None:
            candidate_pairs = self.reverse_candidate_pairs(pair.forward)
        else:
            candidate_pairs = self._replacement_rev_pairs(pair, next_pair)
        for candidate in candidate_pairs:
            if interaction_check(candidate.reverse):
                return candidate
        raise NoSuitablePrimers

    def repair(self, existing_pools: tuple[list[PrimerPair], list[PrimerPair]]) -> None:
        """Repair an existing scheme against a new reference."""

        amplicon_num = 0
        this_existing_pool = existing_pools[self._pool_index]
        other_existing_pool = existing_pools[(self._pool_index + 1) % 2]

        while len(this_existing_pool):
            next_pair: Optional[PrimerPair] = (
                other_existing_pool[0] if len(other_existing_pool) else None
            )
            pair = this_existing_pool.pop(0)
            new_pair = None

            amplicon_num += 1

            fwd_missing = pair.forward.as_kmer() not in self.fwd_kmers
            rev_missing = pair.reverse.as_kmer() not in self.rev_kmers

            if fwd_missing and not rev_missing:
                try:
                    new_pair = self._find_replacement_fwd(pair, this_existing_pool)
                    logger.info("Replaced fwd primer for amplicon {}", amplicon_num)
                except NoSuitablePrimers:
                    pass
            if rev_missing and not fwd_missing:
                try:
                    new_pair = self._find_replacement_rev(
                        pair, next_pair, this_existing_pool
                    )
                    logger.info("Replaced rev primer for amplicon {}", amplicon_num)
                except NoSuitablePrimers:
                    pass

            # Full freedom when replacing both
            if not new_pair and (fwd_missing or rev_missing):
                try:
                    new_pair = self._find_pair(
                        next_pair=next_pair, same_pool_pairs=this_existing_pool
                    )
                    logger.info("Replaced both primers for amplicon {}", amplicon_num)
                except NoSuitablePrimers:
                    logger.warning(
                        "Unable to replace primers for amplicon {}", amplicon_num
                    )

            self._this_pool.append(new_pair if new_pair else pair)
            this_existing_pool = existing_pools[self._pool_index]
            other_existing_pool = existing_pools[(self._pool_index + 1) % 2]
