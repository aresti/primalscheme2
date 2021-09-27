import click
import pathlib
import sys

from typing import Any, TextIO

from .config import Config
from .scheme import DeepScheme
from .utils import check_or_create_outpath, simple_fasta_parser


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
    fasta = list(simple_fasta_parser(input))[0]
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


if __name__ == "__main__":
    main()
