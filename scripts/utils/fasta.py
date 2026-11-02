"""Minimal FASTA reading and writing."""

from typing import Iterator, Tuple


def parse_fasta(handle) -> Iterator[Tuple[str, str]]:
    """Yield ``(record_id, sequence)`` from an open FASTA text handle.

    The id is the first whitespace-delimited token of the header line.
    Raises ``ValueError`` if sequence data appears before any header.
    """
    rec_id, chunks = None, []
    for lineno, line in enumerate(handle, start=1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if rec_id is not None:
                yield rec_id, "".join(chunks)
            header = line[1:].split()
            rec_id, chunks = (header[0] if header else ""), []
        else:
            if rec_id is None:
                raise ValueError(
                    f"Line {lineno}: sequence data before first '>' header")
            chunks.append(line)
    if rec_id is not None:
        yield rec_id, "".join(chunks)


def read_fasta(path) -> list:
    """Read a FASTA file into a list of ``(id, sequence)`` tuples."""
    with open(path) as handle:
        return list(parse_fasta(handle))


def write_fasta(path, records) -> None:
    """Write ``(id, sequence)`` pairs to ``path`` in FASTA format."""
    with open(path, "w") as handle:
        for rec_id, seq in records:
            handle.write(f">{rec_id}\n{seq}\n")
