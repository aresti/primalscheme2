"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the main cli entrypoint for primalscheme.

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
import click
import csv
import hashlib
import os
import pathlib
import sys
import webbrowser

from diskcache import Cache
from loguru import logger
from pathlib import Path
from typing import Any, TextIO

from Bio import SeqIO  # type: ignore

from primalscheme.config import Config
from primalscheme.jackhammer import JackhammerScheme
from primalscheme.overlap import OverlapPriorityScheme
from primalscheme.primer import Kmer, Primer, PrimerPair, digest_to_passing_kmers
from primalscheme.scheme import Scheme


logger = logger.opt(colors=True)


def check_or_create_outpath(path: pathlib.Path, force: bool = False) -> pathlib.Path:
    """
    Check for an existing output dir, require --force to overwrite.
    Create dir if required, return Path obj.
    """
    if path.exists() and not force:
        raise IOError("Directory exists add --force to overwrite")

    path.mkdir(exist_ok=True)
    return path


@click.command()
@click.argument("input", nargs=-1, type=click.File("rt"))
@click.option(
    "--output",
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
    "--min-overlap",
    type=click.IntRange(0),
    default=Config.min_overlap,
)
@click.option("--force/--no-force", default=Config.force)
@click.option("--debug/--no-debug", default=False)
@click.option("-s", "--strategy", type=click.Choice(("o", "j")), default="o")
@click.option(
    "-p",
    "--jackhammer-pools",
    type=click.IntRange(1),
    default=Config.jackhammer_pools,
)
@click.option(
    "-d",
    "--jackhammer-density",
    type=click.FloatRange(0.1, 0.9),
    default=Config.jackhammer_density,
)
@click.option(
    "--prefix",
    type=click.STRING,
    help="Prefix name for your outputs.",
    metavar="<str>",
    default=Config.prefix,
    show_default=True,
)
@click.option(
    "--repair",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=True,
        path_type=pathlib.Path,
    ),
)
@click.option(
    "--repair-interactions/--no-repair-interactions", default=Config.repair_interactions
)
def main(
    input: list[TextIO],
    **kwargs: Any,
) -> None:

    logger.remove()  # Remove default stderr logger
    logger.add(
        sys.stderr,
        colorize=True,
        format="{message}",
        level="DEBUG" if kwargs["debug"] else "INFO",
    )

    if not len(input):
        # noop  - see note in click docs
        # https://click.palletsprojects.com/en/7.x/arguments/?highlight=nargs#variadic-arguments
        sys.exit()

    fastas = [record for file in input for record in SeqIO.parse(file, "fasta")]
    cfg = Config(**kwargs)

    if not cfg.output.is_absolute():
        # Make output path absolute
        cfg.output = Path(os.getcwd()) / cfg.output

    check_or_create_outpath(cfg.output, force=cfg.force)

    logger.add(cfg.output / f"{cfg.prefix}-log.txt", format="{message}", level="DEBUG")

    primary_ref = fastas[0]
    kmers_passing_thermo: list[Kmer] = []

    # Digestion
    cache = Cache("./.primalscheme_cache")
    checksum = hashlib.md5(bytes(primary_ref.seq)).hexdigest()
    if checksum in cache:
        logger.info(
            "Found cached kmers for reference <blue>{ref}</>",
            ref=primary_ref.id,
            min=cfg.primer_size_min,
            max=cfg.primer_size_max,
        )
        kmers_passing_thermo = cache[checksum]
    else:
        logger.info(
            "Digesting reference <blue>{ref}</> into kmers [{min}-{max}]-mers",
            ref=primary_ref.id,
            min=cfg.primer_size_min,
            max=cfg.primer_size_max,
        )

        with click.progressbar(
            length=cfg.primer_size_max - cfg.primer_size_min + 1, label="Digesting"
        ) as pbar:
            kmers_passing_thermo = digest_to_passing_kmers(
                primary_ref.seq, cfg, pbar=pbar
            )
        cache[checksum] = kmers_passing_thermo

    logger.info(
        "Found <blue>{n}</> non-ambiguous kmers passing thermo filter",
        n=len(kmers_passing_thermo),
    )

    if kwargs["strategy"] == "o":
        if cfg.repair:
            existing_pools: tuple[list[PrimerPair], list[PrimerPair]] = ([], [])

            with open(cfg.repair, newline="") as bedfile:
                bedreader = csv.reader(bedfile, delimiter="\t")
                for n, row in enumerate(bedreader):
                    pool = 0 if n % 4 < 2 else 1
                    if n % 2 == 0:
                        fwd = Primer.from_bed_row(row)
                    else:
                        rev = Primer.from_bed_row(row)
                        existing_pools[pool].append(
                            PrimerPair(forward=fwd, reverse=rev)
                        )

            scheme: Scheme = OverlapPriorityScheme(
                primary_ref,
                kmers=kmers_passing_thermo,
                cfg=cfg,
            )
            scheme.repair(existing_pools, fix_interactions=cfg.repair_interactions)
        else:
            with click.progressbar(
                length=len(primary_ref.seq), label="Designing scheme"
            ) as pbar:
                scheme = OverlapPriorityScheme(
                    primary_ref, kmers=kmers_passing_thermo, cfg=cfg, pbar=pbar
                )
                scheme.execute()

    elif kwargs["strategy"] == "j":
        with click.progressbar(
            length=cfg.jackhammer_pools, label="Designing scheme"
        ) as pbar:
            scheme = JackhammerScheme(
                primary_ref, kmers=kmers_passing_thermo, cfg=cfg, pbar=pbar
            )
            scheme.execute()

    # Output
    scheme.write_primer_bed()
    # scheme.write_primer_gff()
    scheme.write_report()

    report_uri = scheme.report_filepath.as_uri()
    logger.success(
        "All done! Scheme created with <blue>{n_pairs}</> amplicons, "
        "<blue>{cov} %</> coverage, <blue>{gaps}</> gaps",
        n_pairs=len(scheme.primer_pairs()),
        cov=scheme.percent_coverage(),
        gaps=scheme.gap_count(),
    )
    logger.info(report_uri)
    webbrowser.open(report_uri)


if __name__ == "__main__":
    main()
