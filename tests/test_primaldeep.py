import random

from primaldeep import __version__
from primaldeep.utils import (
    digest_seq,
    filter_unambiguous_kmers,
    AMBIGUOUS_DNA,
    UNAMBIGUOUS_DNA,
)


def get_random_sequence(length: int = 5000, alphabet: str = UNAMBIGUOUS_DNA) -> str:
    return "".join([random.choice(alphabet) for i in range(length)])


def test_version() -> None:
    assert __version__ == "0.1.0"


def test_digest_seq_output_matches_kmer_size() -> None:
    kmer_size = random.randint(1, 50)
    kmers = digest_seq(get_random_sequence(), kmer_size)
    assert all(len(kmer.seq) == kmer_size for kmer in kmers)


def test_digest_seq_returns_expected_number_of_kmers() -> None:
    kmer_len = random.randint(1, 50)
    seq = get_random_sequence()
    kmers = digest_seq(seq, kmer_len)
    assert len(kmers) == len(seq) - kmer_len + 1


def test_digest_seq_returns_unique_kmers() -> None:
    kmer_size = random.randint(1, 50)
    kmers = digest_seq(get_random_sequence(), kmer_size)
    assert len(set(kmers)) == len(kmers)


def test_filter_unambiguous_kmers() -> None:
    ambig_seq = get_random_sequence(alphabet=AMBIGUOUS_DNA)
    assert not all(b in UNAMBIGUOUS_DNA for b in ambig_seq)

    kmer_size = random.randint(1, 50)
    mixed_kmers = digest_seq(ambig_seq, kmer_size)
    filtered_kmers = list(filter_unambiguous_kmers(mixed_kmers))
    aggregate_filtered_seq = "".join(kmer.seq for kmer in filtered_kmers)

    assert len(mixed_kmers) > len(filtered_kmers)
    assert all(b in UNAMBIGUOUS_DNA for b in aggregate_filtered_seq)
