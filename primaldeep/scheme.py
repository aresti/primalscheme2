import csv
from pathlib import Path
from typing import Sequence

from jinja2 import Environment, PackageLoader, select_autoescape

from primaldeep.config import Config
from primaldeep.datauri import get_data_uri
from primaldeep.dna import reverse_complement, SeqRecordProtocol
from primaldeep.primer import Kmer, Primer, PrimerDirection, PrimerPair


class Scheme:
    def __init__(self, ref: SeqRecordProtocol, kmers: list[Kmer], cfg: Config):
        self.ref = ref
        self.kmers = kmers
        self.cfg = cfg
        self.pools: Sequence[Sequence[PrimerPair]] = []

    def _find_pairs(
        self,
        forward_primer: Primer,
    ) -> Sequence[PrimerPair]:
        """
        For a given forward primer, return all possible PrimerPairs, sorted by reverse
        amplicon length.
        """
        window_left = (
            forward_primer.start + self.cfg.amplicon_size_min - self.cfg.primer_size_max
        )
        window_right = forward_primer.start + self.cfg.amplicon_size_max - 1
        reverse_candidates = [
            k for k in self.kmers if k.start >= window_left and k.end <= window_right
        ]
        candidate_pairs = [
            PrimerPair(
                forward_primer,
                Primer(reverse_complement(r.seq), r.end, PrimerDirection.REVERSE),
            )
            for r in reverse_candidates
        ]
        candidate_pairs = [
            pair for pair in candidate_pairs if len(pair) >= self.cfg.amplicon_size_min
        ]
        candidate_pairs.sort(key=len, reverse=True)
        return candidate_pairs

    @property
    def primer_pairs(self) -> list[PrimerPair]:
        return sorted(
            [pair for pool in self.pools for pair in pool],
            key=lambda p: p.forward.start,
        )

    @property
    def primer_bed_rows(
        self,
    ) -> Sequence[tuple[str, int, int, str, int, str, str]]:
        """Format primer BED rows."""
        rows: list[tuple[str, int, int, str, int, str, str]] = []
        ref_id = self.ref.id
        primer_name_template = self.cfg.prefix + "_{}_{}"

        for n, pair in enumerate(self.primer_pairs):
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
                rows.append(row)

        return rows

    @property
    def primer_bed(self) -> str:
        rows = ["\t".join(map(str, row)) for row in self.primer_bed_rows]
        return "\n".join(rows)

    def write_primer_bed(self) -> None:
        filepath = self.cfg.output / f"{self.cfg.prefix}.primer.bed"

        with filepath.open("w") as fh:
            cw = csv.writer(fh, delimiter="\t")
            cw.writerows(self.primer_bed_rows)

    @property
    def primer_gff(self) -> str:
        rows = []
        ref_id = self.ref.id
        primer_name_template = self.cfg.prefix + "_{}_{}"

        for n, pair in enumerate(self.primer_pairs):
            pool_num = n % len(self.pools) + 1
            color = "#274e13" if pool_num == 1 else "#16537e"
            rows.append(
                [
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
            )
            for direction in PrimerDirection:
                name = primer_name_template.format(n, direction.name)
                p = getattr(pair, direction.name.lower())
                start = p.start if p.direction == PrimerDirection.FORWARD else p.end
                end = p.end if p.direction == PrimerDirection.FORWARD else p.start
                rows.append(
                    [
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
                )
        tsv_rows = ["\t".join(map(str, row)) for row in rows]
        return "##gff-version 3\n" + "\n".join(tsv_rows)

    def write_primer_gff(self) -> None:
        filepath = self.cfg.output / f"{self.cfg.prefix}.primer.gff"

        with filepath.open("w") as fh:
            fh.write(self.primer_gff)

    @property
    def report_filepath(self) -> Path:
        """Filepath for scheme HTML report."""
        return self.cfg.output / f"{self.cfg.prefix}_report.html"

    def write_report(self) -> None:
        """Write HTML report and associated JSON"""

        # Data JSON
        data = {
            "reference": get_data_uri(self.ref.format("fasta")),
            "primerTrack": get_data_uri(self.primer_gff),
        }

        # HTML Template
        jinja_env = Environment(
            loader=PackageLoader("primaldeep"), autoescape=select_autoescape()
        )
        template = jinja_env.get_template("report_template.html")
        rendered = template.render(data=data)

        with self.report_filepath.open("w") as fh:
            fh.write(rendered)
