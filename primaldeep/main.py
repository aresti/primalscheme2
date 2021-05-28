import click
import math

from .parsers import SimpleFastaParser
from typing import Any, TextIO


class Config:
    """A class for managing configuration"""

    amplicon_size = 400
    amplicon_variance = 20
    coverage = 2
    packing = 0.5

    def __init__(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def __str__(self) -> str:
        return str(self.__dict__)


class Region:
    """A defned region within a pool of a scheme"""

    def __init__(
        self,
        ref_slice: str,
        config: Config,
    ):
        self.ref_slice = ref_slice
        self.config = config


class DeepScheme:
    """A scheme with n pools"""

    def __init__(self, fasta: tuple[str, str], config: Config) -> None:
        self.fasta = fasta
        self.config = config
        self.pools: list[list[Region]] = [
            [] for _ in range(math.ceil(config.coverage / config.packing))
        ]
        self.region_size = int(1 / config.packing / config.amplicon_size)
        self.pool_offset = int(self.region_size / len(self.pools))

    @property
    def ref(self) -> str:
        """The reference sequence"""
        return self.fasta[1]

    def __str__(self) -> str:
        return f"DeepScheme: {self.fasta[0]} ({self.config})"


@click.command()
@click.argument("input", type=click.File("rt"))
@click.option(
    "--amplicon-size", type=click.IntRange(100, 2000), default=Config.amplicon_size
)
@click.option(
    "--coverage", type=click.IntRange(0), default=Config.coverage, prompt=True
)
@click.option(
    "--packing", type=click.FloatRange(0.1, 1.0), default=Config.packing, prompt=True
)
def main(
    input: TextIO,
    **kwargs: Any,
) -> None:
    fasta = list(SimpleFastaParser(input))[0]
    config = Config(**kwargs)
    scheme = DeepScheme(fasta, config)
    print(scheme)


if __name__ == "__main__":
    main()
