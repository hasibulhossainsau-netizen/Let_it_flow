"""Pure sequence helpers (no third-party dependencies)."""

import re

_COMPLEMENT = str.maketrans(
    "ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
    "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb",
)


def gc_content(seq: str) -> float:
    """Return the GC fraction (0-1) of ``seq``."""
    if not seq:
        raise ValueError("Cannot compute GC content of an empty sequence")
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def seed_gc(spacer: str, seed_len: int) -> float:
    """GC fraction of the PAM-proximal seed (last ``seed_len`` nt)."""
    return gc_content(spacer[-seed_len:])


def has_homopolymer(seq: str, max_run: int) -> bool:
    """True if ``seq`` contains a single-base run of ``max_run`` or longer."""
    return bool(re.search(r"([ACGT])\1{%d,}" % (max_run - 1), seq.upper()))


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence (IUPAC aware)."""
    return seq.translate(_COMPLEMENT)[::-1]
