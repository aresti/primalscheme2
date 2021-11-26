import re
import parasail  # type: ignore

from primaldeep.dna import CIGAR_REGEX, tm, ThermoConfig


def cig_check(
    regex_matches: list[tuple[str, str]], query: str, cfg: ThermoConfig
) -> bool:
    if regex_matches[0][1] == "=":  # Nothing to extend
        return False

    overlap = ""
    indels = 0
    matches = 0
    mismatches = 0

    for group in regex_matches:  # Iterate over the sequences
        if group[1] == "I" or group[1] == "D":  # Skip I or D ops
            indels += int(group[0])
            continue
        if group[1] == "X":  # Mismatch
            if not overlap:  # Must be 3' anchored
                return False
            mismatches += int(group[0])
        if group[1] == "=":  # Match
            i = sum((indels, matches, mismatches))
            overlap += query[i : i + int(group[0])]
            matches += int(group[0])

        identity = matches / (matches + mismatches)
        ol_tm = tm(overlap, cfg)
        return ol_tm > cfg.dimer_max_tm and identity > cfg.dimer_min_identity

    return False


def parasail_align(seq1: str, seq2: str) -> parasail.Traceback:
    MATRIX = parasail.matrix_create("ACGT", 2, -1)
    OPEN = 10
    EXTEND = 5
    # Semi-Global, do not penalize gaps at beginning and end of both sequences
    trace = parasail.sg_trace_scan(seq1, seq2, OPEN, EXTEND, MATRIX)
    return trace


def interaction_check(seq1: str, seq2: str, cfg: ThermoConfig) -> bool:
    trace = parasail_align(seq1, seq2)
    traceback = trace.get_traceback()
    matches = re.findall(CIGAR_REGEX, trace.cigar.decode.decode())
    return cig_check(matches, traceback.query, cfg) or cig_check(
        matches[::-1], traceback.ref[::-1], cfg
    )
