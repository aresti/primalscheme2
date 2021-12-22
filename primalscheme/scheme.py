"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2021 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the Scheme class.

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

import csv
from operator import attrgetter
from pathlib import Path
from typing import Optional, Protocol, Sequence

from jinja2 import Environment, PackageLoader, select_autoescape

from primalscheme.config import Config
from primalscheme.datauri import get_data_uri
from primalscheme.dna import SeqRecordProtocol
from primalscheme.exceptions import NoReverseWindow
from primalscheme.primer import Kmer, Insert, Primer, PrimerDirection, PrimerPair


class ProgressBar(Protocol):
    def update(self, n_steps: int, current_item: None = None) -> None:
        ...


class Scheme:
    def __init__(
        self,
        ref: SeqRecordProtocol,
        kmers: list[Kmer],
        cfg: Config,
        pbar: ProgressBar = None,
    ):
        self.ref = ref
        self.kmers = kmers
        self.cfg = cfg
        self.pbar = pbar

        self.pools: Sequence[Sequence[PrimerPair]] = []

    def reverse_primer_window_start(
        self, fwd: Primer, next_pair: Optional[PrimerPair] = None
    ) -> int:
        """
        Given a forward Primer, find the start of the window for reverse candidates.
        Min amplicon size and max primer size determine this bound, relative to the
        forward primer.
        """
        if next_pair is None:
            return fwd.start + self.cfg.amplicon_size_min - self.cfg.primer_size_max
        else:
            return next_pair.forward.end

    def reverse_primer_window_end(self, fwd: Primer) -> int:
        """
        Given a forward Primer, find the end of the window for reverse candidates.
        Max amplicon size sets the right bound, relative to the forward primer.
        """
        return fwd.start + self.cfg.amplicon_size_max

    def reverse_primer_window(
        self, fwd: Primer, next_pair: Optional[PrimerPair] = None
    ) -> tuple[int, int]:
        """
        Given a foward Primer, return reference coords describing a window for
        possible reverse candidates.
        """
        start = self.reverse_primer_window_start(fwd, next_pair=next_pair)
        end = self.reverse_primer_window_end(fwd)

        if (end - start) < self.cfg.primer_size_min:
            raise NoReverseWindow
        return (start, end)

    def reverse_candidate_kmers(
        self, fwd: Primer, next_pair: Optional[PrimerPair] = None
    ) -> list[Kmer]:
        """
        Given a forward Primer, return a list of Kmers as candidates for a
        reverse primer, satisfying amplicon size constraints.
        """
        window = self.reverse_primer_window(fwd, next_pair=next_pair)
        return [k for k in self.kmers if k.start >= window[0] and k.end <= window[1]]

    def reverse_candidate_pairs(
        self, fwd: Primer, next_pair: Optional[PrimerPair] = None
    ) -> list[PrimerPair]:
        """
        Given a forward Primer, return a list of candidate PrimerPairs that
        satisfy amplicon size constraints.

        If next_pair is provided, constrain to those pairs that would maintain overlap.
        """
        try:
            reverse_primers = [
                Primer.from_kmer(kmer, PrimerDirection.REVERSE)
                for kmer in self.reverse_candidate_kmers(fwd, next_pair=next_pair)
            ]
        except NoReverseWindow:
            return []

        # Generate all pair combinations
        pairs = [PrimerPair(fwd, rev) for rev in reverse_primers]

        # Remove pairs that generate an amplicon size less than min
        return [pair for pair in pairs if len(pair) >= self.cfg.amplicon_size_min]

    def primer_pairs(self) -> list[PrimerPair]:
        """
        Return a flat list of all PrimerParis in the scheme, sorted by start position.
        """
        return sorted(
            [pair for pool in self.pools for pair in pool],
            key=attrgetter("start"),
        )

    def inserts(self) -> list[Insert]:
        """Return a list of Inserts."""
        return [pair.insert for pair in self.primer_pairs()]

    def _current_amplicon_num(self) -> int:
        return len(self.primer_pairs()) + 1

    def gap_count(self) -> int:
        """The total number of gaps in the scheme."""
        gaps = 0
        inserts = self.inserts()
        last_covered = inserts[0].end
        for insert in inserts[1:]:
            if insert.start > last_covered:
                gaps += 1
            last_covered = insert.end
        return gaps

    def percent_coverage(self) -> float:
        """Coverage %, with respect to the reference."""
        covered_coords = set(
            [x for insert in self.inserts() for x in range(insert.start, insert.end)]
        )
        return round(len(covered_coords) / len(self.ref.seq) * 100, 2)

    def primer_bed_rows(
        self,
    ) -> Sequence[tuple[str, ...]]:
        """
        Return BED formatted rows for all primers in the scheme.
        """
        rows = []
        ref_id = self.ref.id
        primer_name_template = self.cfg.prefix + "_{}_{}"

        for n, pair in enumerate(self.primer_pairs()):
            pool_num = n % len(self.pools) + 1
            for direction in PrimerDirection:
                name = primer_name_template.format(
                    n + 1, "LEFT" if direction == PrimerDirection.FORWARD else "RIGHT"
                )
                p: Primer = getattr(pair, direction.name.lower())
                start = p.start
                end = p.end
                row = (
                    ref_id,
                    start,
                    end,
                    name,
                    pool_num,
                    direction.value,
                    p.seq,
                )
                rows.append(tuple(map(str, row)))

        return rows

    def write_primer_bed(
        self,
    ) -> None:
        """
        Write a BED file describing all primers in the scheme.
        """
        filepath = self.cfg.output / f"{self.cfg.prefix}.primer.bed"

        with filepath.open("w") as fh:
            cw = csv.writer(fh, delimiter="\t")
            cw.writerows(self.primer_bed_rows())

    def amplicon_primer_gff_rows(self) -> list[tuple[str, ...]]:

        rows = []
        ref_id = self.ref.id
        primer_name_template = self.cfg.prefix + "_{}_{}"

        # Amplicons
        for n, pair in enumerate(self.primer_pairs()):
            pool_num = n % len(self.pools) + 1
            color = "#274e13" if pool_num == 1 else "#16537e"
            row = [
                ref_id,
                "primalscheme",
                "mRNA",
                pair.start + 1,
                pair.end,  # GFF is closed
                ".",
                ".",
                ".",
                f"ID=AMP{n+1};Color={color}",
            ]
            rows.append(tuple(map(str, row)))

            # Primers
            for direction in PrimerDirection:
                name = primer_name_template.format(n, direction.name)
                p = getattr(pair, direction.name.lower())
                start = p.start
                end = p.end
                row = [
                    ref_id,
                    "primalscheme",
                    "exon",
                    start + 1,
                    end,
                    ".",
                    ".",
                    ".",
                    f"ID={name};Name={name};Parent=AMP{n+1};Note=POOL_{pool_num}",
                ]
                rows.append(tuple(map(str, row)))

        return rows

    def amplicon_primer_gff(self) -> str:
        """
        Return GFF3 formatted multiline str for all amplicons and associated primers in
        the scheme.
        """
        tsv_rows = ["\t".join(map(str, row)) for row in self.amplicon_primer_gff_rows()]
        return "##gff-version 3\n" + "\n".join(tsv_rows)

    def write_primer_gff(self) -> None:
        """
        Write a GFF3 file for all amplicons and associated primers in the scheme.
        """
        filepath = self.cfg.output / f"{self.cfg.prefix}.primer.gff"

        with filepath.open("w") as fh:
            fh.write(self.amplicon_primer_gff())

    @property
    def report_filepath(self) -> Path:
        """Filepath for scheme HTML report."""
        return self.cfg.output / f"{self.cfg.prefix}_report.html"

    def write_report(self) -> None:
        """
        Write HTML report and associated JSON
        """

        # Data JSON
        data = {
            "reference": get_data_uri(self.ref.format("fasta")),
            "primerTrack": get_data_uri(self.amplicon_primer_gff()),
        }

        # HTML Template
        jinja_env = Environment(
            loader=PackageLoader("primalscheme"), autoescape=select_autoescape()
        )
        template = jinja_env.get_template("report_template.html")
        rendered = template.render(data=data)

        with self.report_filepath.open("w") as fh:
            fh.write(rendered)
