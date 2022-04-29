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
        """Ideal spacing between amplicons, for the given density."""
        return int(self.cfg.amplicon_size_target / self.cfg.jackhammer_density)

    def execute(self) -> None:
        """Execute the scheme design."""
        for n, pool in enumerate(self.pools):
            offset = int(n * self.spacing / len(self.pools))
            pool.extend(JackhammerPool(scheme=self, num=n, offset=offset).execute())

            if self.pbar:
                self.pbar.update(1)


class JackhammerPool:
    def __init__(self, scheme: JackhammerScheme, num: int, offset: int):
        self.scheme = scheme
        self.num = num
        self.offset = offset

        self.pairs: list[PrimerPair] = []
        self._amplicon_num = 0

    @property
    def _prev_pair(self) -> Optional[PrimerPair]:
        """Previous pair in the pool."""
        return self.pairs[-1] if len(self.pairs) else None

    @property
    def primers(self) -> list[Primer]:
        """A flat list of the primers in this pool."""
        return [p for pp in self.pairs for p in (pp.forward, pp.reverse)]

    @property
    def _next_optimal_start_position(self) -> int:
        """Next optimal start position to maintain ideal spacing."""
        return self.offset + self._amplicon_num * self.scheme.spacing

    @property
    def _ref_end(self) -> int:
        return len(self.scheme.ref.seq)

    @property
    def _last_start_for_end(self) -> int:
        """
        The last start position of a forward primer, forming a potential amplicon of
        min size, that would not extend past the end of the reference.
        """
        return self._ref_end - self.scheme.cfg.amplicon_size_min

    @property
    def _last_start_for_cell(self) -> int:
        """
        The last start position of a forward primer, forming a potential amplicon of
        min size, that would not breach the next optimum start position.
        """
        return (
            (self._amplicon_num + 1) * self.scheme.spacing
            - self.scheme.cfg.amplicon_size_min
            + self.offset
        )

    @property
    def _last_useful_start(self) -> int:
        """
        The last useful start position of a forward primer for the current amplicon.
        """
        return min(self._last_start_for_cell, self._last_start_for_end)

    def _available_fwd_kmers(self) -> list[Kmer]:
        """All available (non-crashing) forward kmers."""
        if self._prev_pair is None:
            # Empty scheme
            return self.scheme.kmers

        start_from = self._prev_pair.reverse.end

        return [
            k
            for k in self.scheme.kmers
            if k.start >= start_from and k.start <= self._last_useful_start
        ]

    def _sorted_fwd_kmers(self) -> list[Kmer]:
        """Forward kmers, sorted by deviation from optimal start position."""
        available_kmers = self._available_fwd_kmers()
        return sorted(
            available_kmers,
            key=lambda k: abs(self._next_optimal_start_position - k.start),
        )

    def _interaction_checker_factory(
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

        def inner_func(primer: Primer, fwd: Primer = None) -> bool:
            # Hairpin
            if primer.forms_hairpin(self.scheme.cfg):
                logger.debug(f"Primer {primer}: hairpin predicted.")
                return False

            # Self
            if primer.interacts_with([primer], self.scheme.cfg, verbose=verbose):
                logger.debug(f"Primer {primer}: interacts with self.")
                return False

            # Previous pairs
            if primer.interacts_with(check_primers, self.scheme.cfg, verbose=verbose):
                logger.debug(
                    f"Primer {primer}: interacts with another primer in the pool"
                )
                return False

            # Forward primer
            if fwd and primer.interacts_with([fwd], self.scheme.cfg, verbose=verbose):
                logger.debug(
                    f"Reverse primer {primer}: interacts with the selected forward primer."
                )
                return False

            return True

        return inner_func

    def _sorted_reverse_candidate_pairs(self, fwd: Primer) -> list[PrimerPair]:
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
        passes_interaction_checks = self._interaction_checker_factory(
            existing_pairs=same_pool_pairs
        )

        for fwd_kmer in self._sorted_fwd_kmers():
            fwd_primer = Primer.from_kmer(fwd_kmer, PrimerDirection.FORWARD)

            # Find pairs
            try:
                candidate_pairs = self._sorted_reverse_candidate_pairs(fwd_primer)
            except NoSuitablePrimers:
                continue

            # Check interactions
            if not passes_interaction_checks(fwd_primer):
                continue
            for pair in candidate_pairs:
                if passes_interaction_checks(pair.reverse, fwd=fwd_primer):
                    return pair

        raise NoSuitablePrimers

    def execute(self) -> list[PrimerPair]:
        """Execute the design of the pool."""
        while True:
            try:
                pair = self._find_pair()
                pair.pool = self.num
                self.pairs.append(pair)
            except NoSuitablePrimers:
                if self._last_useful_start == self._last_start_for_end:
                    # End of scheme
                    break
            self._amplicon_num += 1

        return self.pairs
