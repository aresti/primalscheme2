"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2021 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the JackhammerScheme class.

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
from typing import Callable, Optional

from Bio import SeqRecord  # type: ignore

from primalscheme.exceptions import NoSuitablePrimers
from primalscheme.primer import Kmer, Primer, PrimerDirection, PrimerPair
from primalscheme.scheme import Scheme, ProgressBar
from primalscheme.config import Config

logger = logger.opt(colors=True)


class JackhammerScheme(Scheme):
    def __init__(
        self,
        ref: SeqRecord,
        kmers: list[Kmer],
        cfg: Config,
        pbar: ProgressBar = None,
    ):
        super().__init__(ref, kmers, cfg, pbar=pbar)

        self.pools: tuple[list[PrimerPair], ...] = tuple(
            [] for _ in range(self.cfg.jackhammer_pools)
        )

    @property
    def spacing(self) -> int:
        return int(self.cfg.amplicon_size_target / self.cfg.jackhammer_density)

    def execute(self) -> None:
        for n, pool in enumerate(self.pools):
            offset = int(n * self.spacing / len(self.pools))
            pool.extend(JackhammerPool(scheme=self, offset=offset).execute())

            if self.pbar:
                self.pbar.update(1)


class JackhammerPool:
    def __init__(self, scheme: JackhammerScheme, offset: int):
        self.scheme = scheme
        self.offset = offset

        self.pairs: list[PrimerPair] = []

    @property
    def _prev_pair(self) -> Optional[PrimerPair]:
        """
        Return the previous pair.
        """
        return self.pairs[-1] if len(self.pairs) else None

    @property
    def _amplicon_num(self) -> int:
        return len(self.pairs)

    @property
    def primers(self) -> list[Primer]:
        """
        Return a flat list of the primers in this pool.
        """
        return [p for pp in self.pairs for p in (pp.forward, pp.reverse)]

    @property
    def _next_optimal_start_position(self) -> int:
        return self.offset + self._amplicon_num * self.scheme.spacing

    def _available_fwd_kmers(self) -> list[Kmer]:
        """Return all available (non-crashing) forward kmers."""
        if self._prev_pair is None:
            # Empty scheme
            return self.scheme.kmers

        start_from = self._prev_pair.reverse.end

        # We don't want fwd kmers that would result in the amplicon extending beyond
        # the next optimum start position
        last_start_for_cell = (
            (self._amplicon_num + 1) * self.scheme.spacing
            - self.scheme.cfg.amplicon_size_min
            + self.offset
        )

        # Also, we don't want to fwd kmers that can't possibly make a large enough
        # amplicon
        ref_end = len(self.scheme.ref.seq)
        last_start_for_end = ref_end - self.scheme.cfg.amplicon_size_min

        last_useful_start = min(last_start_for_cell, last_start_for_end)

        return [
            k
            for k in self.scheme.kmers
            if k.start >= start_from and k.start <= last_useful_start
        ]

    def _sorted_fwd_kmers(self) -> list[Kmer]:
        available_kmers = self._available_fwd_kmers()
        return sorted(
            available_kmers,
            key=lambda k: abs(self._next_optimal_start_position - k.start),
        )

    def interaction_checker_factory(
        self, existing_pairs: Optional[list[PrimerPair]] = None, verbose: bool = False
    ) -> Callable[[Primer], bool]:
        """
        Return a function that performs the required interaction checks for a candidate
        Primer against all previously selected primers in the current pool.
        Returns True where no interactions are found.
        """

        check_primers = self.primers
        if existing_pairs is not None:
            check_primers += [
                p for pp in existing_pairs for p in (pp.forward, pp.reverse)
            ]

        def inner_func(primer: Primer) -> bool:
            if primer.forms_hairpin(self.scheme.cfg):
                logger.trace(f"Primer {primer}: hairpin predicted.")
                return False
            return not (
                primer.interacts_with([primer], self.scheme.cfg, verbose=verbose)
                or primer.interacts_with(
                    check_primers, self.scheme.cfg, verbose=verbose
                )
            )

        return inner_func

    def sorted_reverse_candidate_pairs(self, fwd: Primer) -> list[PrimerPair]:
        """
        Given a forward Primer, return a list of candidate PrimerPairs that
        satisfy amplicon size constraints, sorted by amplicon size deviation from mean.
        """
        pairs = self.scheme.reverse_candidate_pairs(fwd)

        if not len(pairs):
            raise NoSuitablePrimers

        # Sort by deviation from target
        pairs.sort(
            key=lambda p: abs(self.scheme.cfg.amplicon_size_target - p.amplicon_size)
        )
        return pairs

    def _find_pair(
        self, same_pool_pairs: Optional[list[PrimerPair]] = None
    ) -> PrimerPair:
        """Find the next PrimerPair for the scheme."""
        passes_interaction_checks = self.interaction_checker_factory(
            existing_pairs=same_pool_pairs
        )

        for fwd_kmer in self._sorted_fwd_kmers():
            fwd_primer = Primer.from_kmer(fwd_kmer, PrimerDirection.FORWARD)

            # Find pairs
            try:
                candidate_pairs = self.sorted_reverse_candidate_pairs(fwd_primer)
            except NoSuitablePrimers:
                continue

            # Check interactions
            if not passes_interaction_checks(fwd_primer):
                continue
            for pair in candidate_pairs:
                if passes_interaction_checks(pair.reverse):
                    return pair

        raise NoSuitablePrimers

    def execute(self) -> list[PrimerPair]:
        while True:
            try:
                self.pairs.append(self._find_pair())
            except NoSuitablePrimers:
                break

        return self.pairs
