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
import os
import pathlib
import sys
import webbrowser

from typing import Any, TextIO
from operator import attrgetter
from pathlib import Path

from Bio import SeqIO  # type: ignore

from primaldeep.config import Config
from primaldeep.primer import (
    Kmer,
    digest_seq,
    kmer_thermo_check,
    filter_unambiguous_kmers,
)
from primaldeep.overlap import OverlapPriorityScheme


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
def main(
    input: list[TextIO],
    **kwargs: Any,
) -> None:

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
    passing_kmers: list[Kmer] = []

    # Digestion
    kmer_sizes = range(cfg.primer_size_min, cfg.primer_size_max + 1)

    with click.progressbar(kmer_sizes, label="Digesting reference") as sizes:
        for size in sizes:
            digested = [
                Kmer(k.seq, k.start) for k in digest_seq(str(primary_ref.seq), size)
            ]
            unambiguous = filter_unambiguous_kmers(digested)
            passing_kmers.extend(k for k in unambiguous if kmer_thermo_check(k, cfg))

    passing_kmers.sort(key=attrgetter("start"))
    click.echo(f"Found {len(passing_kmers)} non-ambiguous kmers passing thermo filter")

    if kwargs["strategy"] == "o":
        with click.progressbar(
            label="Designing scheme", length=len(primary_ref.seq)
        ) as pbar:
            scheme = OverlapPriorityScheme(
                primary_ref,
                kmers=passing_kmers,
                cfg=cfg,
                pbar=pbar,
            )

        scheme.write_primer_bed()
        scheme.write_primer_gff()
        scheme.write_report()

        report_uri = scheme.report_filepath.as_uri()
        click.echo(report_uri)
        webbrowser.open(report_uri)


if __name__ == "__main__":
    main()
