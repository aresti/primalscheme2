from __future__ import annotations

import pathlib

from primer3 import (
    calcTm as p3_calcTm,
    calcHairpin as p3_calcHairpin,
)
from typing import Iterable, TextIO, TYPE_CHECKING

from .config import Config

if TYPE_CHECKING:
    from .types import Primer


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


def primer_thermo_filter(primer: "Primer", config: Config) -> bool:
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


# From BioPython
def simple_fasta_parser(handle: TextIO) -> Iterable[tuple[str, str]]:
    """Iterate over Fasta records as string tuples.
    Arguments:
     - handle - input stream opened in text mode
    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in simple_fasta_parser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines: list[str] = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")
