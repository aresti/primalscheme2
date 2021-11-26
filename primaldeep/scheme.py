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
from pathlib import Path
from typing import Sequence, Protocol

from jinja2 import Environment, PackageLoader, select_autoescape

from primaldeep.config import Config
from primaldeep.datauri import get_data_uri
from primaldeep.dna import reverse_complement, SeqRecordProtocol
from primaldeep.primer import Kmer, Primer, PrimerDirection, PrimerPair


class ProgressBar(Protocol):
    def update(self, n_steps: int, current_item: None = None) -> None:
        ...


class Scheme:
    def __init__(
        self, ref: SeqRecordProtocol, kmers: list[Kmer], cfg: Config, pbar: ProgressBar
    ):
        self.ref = ref
        self.kmers = kmers
        self.cfg = cfg
        self.pbar = pbar

        self.pools: Sequence[Sequence[PrimerPair]] = []

    def reverse_primer_window_start(self, fwd: Primer) -> int:
        """
        Given a forward Primer, find the start of the window for reverse candidates.
        Min amplicon size and max primer size determine this bound, relative to the
        forward primer.
        """
        return fwd.start + self.cfg.amplicon_size_min - self.cfg.primer_size_max

    def reverse_primer_window_end(self, fwd: Primer) -> int:
        """
        Given a forward Primer, find the end of the window for reverse candidates.
        Max amplicon size sets the right bound, relative to the forward primer.
        """
        return fwd.start + self.cfg.amplicon_size_max - 1

    def reverse_primer_window(self, fwd: Primer) -> tuple[int, int]:
        """
        Given a foward Primer, return reference coords describing a window for
        possible reverse candidates.
        """
        return (
            self.reverse_primer_window_start(fwd),
            self.reverse_primer_window_end(fwd),
        )

    def reverse_candidate_kmers(self, fwd: Primer) -> list[Kmer]:
        """
        Given a forward Primer, return a list of Kmers as candidates for a
        reverse primer, satisfying amplicon size constraints.
        """
        window = self.reverse_primer_window(fwd)
        return [k for k in self.kmers if k.start >= window[0] and k.end <= window[1]]

    def reverse_candidate_pairs(self, fwd: Primer) -> list[PrimerPair]:
        """
        Given a forward Primer, return a list of candidate PrimerPairs that
        satisfy amplicon size constraints, sorted highest-length first.
        """
        pairs = [
            PrimerPair(
                fwd,
                Primer(reverse_complement(kmer.seq), kmer.end, PrimerDirection.REVERSE),
            )
            for kmer in self.reverse_candidate_kmers(fwd)
        ]

        # Remove pairs that generate an amplicon size less than min
        pairs = [pair for pair in pairs if len(pair) >= self.cfg.amplicon_size_min]

        # Reverse length sort
        pairs.sort(key=len, reverse=True)
        return pairs

    def primer_pairs(self) -> list[PrimerPair]:
        """
        Return a flat list of all PrimerParis in the scheme, sorted by start position.
        """
        return sorted(
            [pair for pool in self.pools for pair in pool],
            key=lambda p: p.forward.start,
        )

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
                name = primer_name_template.format(n, direction.name)
                p: Primer = getattr(pair, direction.name.lower())
                start = p.start if p.direction == PrimerDirection.FORWARD else p.end
                end = p.end if p.direction == PrimerDirection.FORWARD else p.start
                row = (
                    ref_id,
                    start,
                    end + 1,  # BED is half-open
                    name,
                    pool_num,
                    direction.value,
                    p.seq,
                )
                rows.append(tuple(map(str, row)))

        return rows

    def write_primer_bed(self) -> None:
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
                pair.end + 1,
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
                start = p.start if p.direction == PrimerDirection.FORWARD else p.end
                end = p.end if p.direction == PrimerDirection.FORWARD else p.start
                row = [
                    ref_id,
                    "primalscheme",
                    "exon",
                    start + 1,
                    end + 1,
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
            loader=PackageLoader("primaldeep"), autoescape=select_autoescape()
        )
        template = jinja_env.get_template("report_template.html")
        rendered = template.render(data=data)

        with self.report_filepath.open("w") as fh:
            fh.write(rendered)
