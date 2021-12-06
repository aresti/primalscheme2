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
import os
import pathlib
import sys
import webbrowser

from loguru import logger
from operator import attrgetter
from pathlib import Path
from typing import Any, TextIO

from Bio import SeqIO  # type: ignore

from primaldeep.config import Config
from primaldeep.overlap import OverlapPriorityScheme
from primaldeep.primer import (
    Kmer,
    Primer,
    PrimerPair,
    digest_seq,
    filter_unambiguous_kmers,
)


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
@click.option("--strategy", type=click.Choice(("o", "p")), default="o")
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
def main(
    input: list[TextIO],
    **kwargs: Any,
) -> None:

    logger.remove()  # Remove default stderr logger
    logger.add(sys.stderr, colorize=True, format="{message}")

    if not len(input):
        # noop  - see note in click docs
        # https://click.palletsprojects.com/en/7.x/arguments/?highlight=nargs#variadic-arguments
        sys.exit()

    fastas = [record for file in input for record in SeqIO.parse(file, "fasta")]
    cfg = Config(**kwargs)

    if not cfg.output.is_absolute():
        # Make output path absolute
        cfg.output = Path(os.getcwd()) / cfg.output

    primary_ref = fastas[0]
    kmers_passing_thermo: list[Kmer] = []
    fwd_kmers: list[Kmer] = []
    rev_kmers: list[Kmer] = []

    # Digestion
    kmer_sizes = range(cfg.primer_size_min, cfg.primer_size_max + 1)
    logger.info(
        "Digesting reference <blue>{ref}</> into kmers [{min}-{max}]-mers",
        ref=primary_ref.id,
        min=cfg.primer_size_min,
        max=cfg.primer_size_max,
    )

    with click.progressbar(
        kmer_sizes, file=sys.stderr, label=f"{cfg.primer_size_min}-mers"
    ) as kmer_sizes_bar:
        for size in kmer_sizes_bar:
            kmer_sizes_bar.label = f"{size}-mers"
            digested = [
                Kmer(k.seq, k.start) for k in digest_seq(str(primary_ref.seq), size)
            ]
            unambiguous = filter_unambiguous_kmers(digested)
            kmers_passing_thermo.extend(
                k for k in unambiguous if k.passes_thermo_checks(cfg)
            )

    kmers_passing_thermo.sort(key=attrgetter("start"))

    logger.info(
        "Found <blue>{n}</> non-ambiguous kmers passing thermo filter",
        n=len(kmers_passing_thermo),
    )

    # Hairpin check
    with click.progressbar(
        kmers_passing_thermo, file=sys.stderr, label="Screening hairpin interactions"
    ) as hairpin_bar:
        for fwd_kmer in hairpin_bar:
            if fwd_kmer.passes_hairpin_check(cfg):
                fwd_kmers.append(fwd_kmer)
            rev_kmer = fwd_kmer.as_reverse()
            if rev_kmer.passes_hairpin_check(cfg):
                rev_kmers.append(rev_kmer)

    logger.info(
        "Found <blue>{n_fwd}</> fwd, <blue>{n_rev}</> rev kmers passing hairpin check",
        n_fwd=len(fwd_kmers),
        n_rev=len(rev_kmers),
    )

    if kwargs["strategy"] == "o":
        scheme = OverlapPriorityScheme(
            primary_ref,
            fwd_kmers=fwd_kmers,
            rev_kmers=rev_kmers,
            cfg=cfg,
        )

        if cfg.repair is None:
            scheme.execute()

        else:
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

            scheme.repair(existing_pools)

        scheme.write_primer_bed()
        scheme.write_primer_gff()
        scheme.write_report()

        report_uri = scheme.report_filepath.as_uri()
        logger.info(report_uri)
        webbrowser.open(report_uri)


if __name__ == "__main__":
    main()
