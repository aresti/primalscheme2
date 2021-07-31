import click
import math
import pathlib
import sys
import csv

from collections import namedtuple
from enum import Enum
from .parsers import SimpleFastaParser
from typing import Any, Iterable, Optional, TextIO

from primer3 import (
    calcTm as p3_calcTm,
    calcHairpin as p3_calcHairpin,
)


Feature = namedtuple("Feature", "start end")  # 0-based, always forward, closed
Kmer = namedtuple("Kmer", "seq start")

UNAMBIGUOUS_DNA = "ACGT"
AMBIGUOUS_DNA = UNAMBIGUOUS_DNA + "RYWSMKHBVDN-"


class NoSuitablePrimersError(Exception):
    """Unable to find any passing primers within the region."""

    pass


def digest_seq(seq: str, kmer_size: int) -> list[Kmer]:
    """Digest a sequence into kmers."""
    return [Kmer(seq[i : i + kmer_size], i) for i in range((len(seq) - kmer_size) + 1)]


def filter_unambiguous_kmers(kmers: Iterable[Kmer]) -> Iterable[Kmer]:
    """Filter out kmers with ambiguity codes in their seq"""
    return filter(lambda k: all(b.upper() in UNAMBIGUOUS_DNA for b in k.seq), kmers)


class Config:
    """
    PrimalDeep configuration.
    Class properties are defaults, can be overriden
    on instantiation (and will shadow class defaults)
    """

    output = pathlib.Path("./output")
    prefix = "scheme"
    force = False

    amplicon_size_min = 380
    amplicon_size_max = 420
    coverage = 2
    packing = 0.5
    high_gc = False

    primer_size_default_min = 19
    primer_size_default_max = 34
    primer_size_default_opt = 22
    primer_size_hgc_min = 17
    primer_size_hgc_max = 30
    primer_size_hgc_opt = 20

    primer_gc_default_min = 30
    primer_gc_default_max = 55
    primer_gc_hgc_min = 40
    primer_gc_hgc_max = 65
    primer_tm_min = 59.5
    primer_tm_max = 62.5
    primer_hairpin_th_max = 47.0
    primer_homopolymer_max = 5

    mv_conc = 100.0
    dv_conc = 2.0
    dntp_conc = 0.8
    dna_conc = 15.0

    def __init__(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

        if self.high_gc:
            self.primer_size_min = self.primer_size_hgc_min
            self.primer_size_max = self.primer_size_hgc_max
            self.primer_size_opt = self.primer_size_hgc_opt
            self.primer_gc_min = self.primer_gc_hgc_min
            self.primer_gc_max = self.primer_gc_hgc_max
        else:
            self.primer_size_min = self.primer_size_default_min
            self.primer_size_max = self.primer_size_default_max
            self.primer_size_opt = self.primer_size_default_opt
            self.primer_gc_min = self.primer_gc_default_min
            self.primer_gc_max = self.primer_gc_default_max

        # Derived
        self.amplicon_size_var = math.ceil(
            (self.amplicon_size_max - self.amplicon_size_min) / 2
        )
        self.amplicon_size_avg = self.amplicon_size_max - self.amplicon_size_var
        self.insert_size_avg = self.amplicon_size_avg - (2 * self.primer_size_opt)
        self.insert_size_max = self.amplicon_size_max - (2 * self.primer_size_min)

    def __str__(self) -> str:
        return str(self.__dict__)


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


class DeepScheme:
    """A scheme with n pools"""

    def __init__(self, fasta: tuple[str, str], cfg: Config) -> None:
        self.fasta = fasta
        self.cfg = cfg

        self.ref = fasta[1]
        self.ref_id = fasta[0]
        max_insert_percent = cfg.amplicon_size_max / cfg.insert_size_max
        self.pools: list[list[Region]] = [
            []
            for _ in range(math.ceil((cfg.coverage / cfg.packing) * max_insert_percent))
        ]
        self.region_size = int(1 / cfg.packing * cfg.amplicon_size_max)
        self.pool_offset = int(self.region_size / len(self.pools))

    def run(self) -> None:
        for n, pool in enumerate(self.pools):
            offset = n * self.pool_offset
            print(f"Pool {n}, offset {offset}")
            for i, start in enumerate(range(offset, len(self.ref), self.region_size)):
                print(f"Region {i} starting at {start}")
                ref_slice = self.ref[
                    start : min((start + self.region_size), len(self.ref))
                ]
                if len(ref_slice) >= self.cfg.amplicon_size_max * 2:
                    region = Region(
                        ref_slice,
                        start,
                        self.cfg,
                    )
                    print(region)
                    pool.append(region)
                else:
                    print("Pool end")

    def __str__(self) -> str:
        return (
            f"DeepScheme: {self.fasta[0]}, region size: {self.region_size} ({self.cfg})"
        )

    def write_primer_bed(self) -> None:
        """Write primer BED file."""
        filepath = self.cfg.output / f"{self.cfg.prefix}.primer.bed"

        rows = []
        ref_id = self.ref_id
        primer_name_template = "POOL_{}_REGION_{}_{}"

        for pool_num, pool in enumerate(self.pools):
            for region_num, region in enumerate(pool):
                for direction in PrimerDirection:
                    name = primer_name_template.format(
                        pool_num, region_num, direction.name
                    )
                    p = getattr(region.top_pair, direction.name.lower())
                    start = p.start if p.direction == PrimerDirection.FORWARD else p.end
                    end = p.end if p.direction == PrimerDirection.FORWARD else p.start
                    rows.append(
                        [
                            ref_id,
                            start,
                            end + 1,  # BED is half-open
                            name,
                            pool_num,
                            direction.value,
                            p.seq,
                        ]
                    )

            with filepath.open("w") as fh:
                cw = csv.writer(fh, delimiter="\t")
                cw.writerows(rows)

    def write_insert_bed(self) -> None:
        """Write insert BED file."""
        filepath = self.cfg.output / f"{self.cfg.prefix}.insert.bed"

        rows = []
        ref_id = self.ref_id
        insert_name_template = "POOL_{}_REGION_{}"

        for pool_num, pool in enumerate(self.pools):
            for region_num, region in enumerate(pool):
                name = insert_name_template.format(pool_num, region_num)
                insert = region.insert
                rows.append(
                    [
                        ref_id,
                        insert.start,
                        insert.end + 1,  # BED is half-open
                        name,
                        "+",
                    ]
                )

            with filepath.open("w") as fh:
                cw = csv.writer(fh, delimiter="\t")
                cw.writerows(rows)


@click.command()
@click.argument("input", type=click.File("rt"))
@click.argument(
    "output",
    type=click.Path(
        exists=False,
        file_okay=False,
        dir_okay=True,
        writable=True,
        path_type=pathlib.Path,
    ),
    default=Config.output,
)
@click.option(
    "--amplicon-size-min",
    type=click.IntRange(100, 2000),
    default=Config.amplicon_size_min,
)
@click.option(
    "--amplicon-size-max",
    type=click.IntRange(100, 2000),
    default=Config.amplicon_size_max,
)
@click.option(
    "--coverage", type=click.IntRange(0), default=Config.coverage, prompt=True
)
@click.option("--force/--no-force", default=Config.force)
@click.option(
    "--packing", type=click.FloatRange(0.1, 1.0), default=Config.packing, prompt=True
)
@click.option(
    "--prefix",
    type=click.STRING,
    help="Prefix name for your outputs.",
    metavar="<str>",
    default=Config.prefix,
    show_default=True,
)
def main(
    input: TextIO,
    **kwargs: Any,
) -> None:
    fasta = list(SimpleFastaParser(input))[0]
    cfg = Config(**kwargs)
    scheme = DeepScheme(fasta, cfg)
    print(scheme)

    # Validate output path
    try:
        cfg.output = check_or_create_outpath(cfg.output, force=cfg.force)
    except IOError as e:
        click.echo(
            click.style(
                f"Error: {e}",
                fg="red",
            )
        )
        sys.exit(1)

    scheme.run()
    scheme.write_primer_bed()
    scheme.write_insert_bed()


def calc_gc(seq: str) -> float:
    """Calculate percent GC for a sequence."""
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def calc_tm(seq: str, config: Config) -> float:
    """Calculate Tm for a sequence."""
    return p3_calcTm(
        seq,
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dntp_conc=config.dntp_conc,
        dna_conc=config.dna_conc,
    )


def calc_max_homo(seq: str) -> int:
    """Calculate max homopolymer length for a sequence."""
    running_homo = 0
    max_homo = 0
    prev_base = None
    for base in seq:
        if base == prev_base:
            running_homo += 1
        else:
            prev_base = base
            running_homo = 1
        if running_homo > max_homo:
            max_homo = running_homo
    return max_homo


def calc_hairpin(seq: str, config: Config) -> object:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Return primer3 thermo object.
    """
    return p3_calcHairpin(
        seq,
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dntp_conc=config.dntp_conc,
        dna_conc=config.dna_conc,
    )


def primer_thermo_filter(primer: Primer, config: Config) -> bool:
    """Hard filter for candidate primers."""
    return (
        (config.primer_gc_min <= primer.gc <= config.primer_gc_max)
        and (config.primer_tm_min <= primer.tm <= config.primer_tm_max)
        and (primer.hairpin <= config.primer_hairpin_th_max)
        and (primer.max_homo <= config.primer_homopolymer_max)
    )


def reverse(seq: str) -> str:
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq: str) -> str:
    """Returns the complement of seq"""
    pairs = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "N": "N",
    }  # todo IUPAC fallback to same
    return "".join(pairs[base] for base in seq)


def reverse_complement(seq: str) -> str:
    """ "Returns the reverse complement of seq"""
    return complement(reverse(seq))


def check_or_create_outpath(path: pathlib.Path, force: bool = False) -> pathlib.Path:
    """
    Check for an existing output dir, require --force to overwrite.
    Create dir if required, return Path obj.
    """
    if path.exists() and not force:
        raise IOError("Directory exists add --force to overwrite")

    path.mkdir(exist_ok=True)
    return path


if __name__ == "__main__":
    main()
