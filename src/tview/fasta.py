"""FASTA parsing and panel construction.

Parses simple FASTA files and builds alignment ``Panel`` objects where the
first sequence is treated as the reference row.

Examples:
    >>> from pathlib import Path
    >>> from tview.fasta import read_fasta, fasta_panel
"""

from __future__ import annotations

from pathlib import Path

from tview.models import Panel


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Parse a FASTA file into a list of (name, sequence) tuples.

    Args:
        path: Path to the FASTA file.

    Returns:
        List of (header_name, concatenated_sequence) tuples.

    Examples:
        >>> import tempfile; from pathlib import Path
        >>> d = Path(tempfile.mkdtemp())
        >>> fasta = d / "test.fa"
        >>> _ = fasta.write_text(">seq1\\nACGT\\n>seq2\\nTGCA\\n")
        >>> read_fasta(fasta)
        [('seq1', 'ACGT'), ('seq2', 'TGCA')]
        >>> read_fasta(d / "empty.fa")
        Traceback (most recent call last):
            ...
        FileNotFoundError: ...
    """
    seqs: list[tuple[str, str]] = []
    name: str | None = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def fasta_panel(
    path: str | Path,
    col_start: int | None = None,
    col_end: int | None = None,
) -> Panel:
    """Build a Panel from an aligned FASTA where the first sequence is the reference.

    Args:
        path: Path to the aligned FASTA file.
        col_start: 1-based inclusive start column for slicing the alignment.
        col_end: 1-based inclusive end column for slicing the alignment.

    Returns:
        A Panel with reference row, sequence rows, and column labels.

    Raises:
        ValueError: If the FASTA file contains no sequences.

    Examples:
        >>> import tempfile; from pathlib import Path
        >>> d = Path(tempfile.mkdtemp())
        >>> fasta = d / "aln.fa"
        >>> _ = fasta.write_text(">ref\\nACGT\\n>read1\\nACTT\\n>read2\\nA-GT\\n")
        >>> p = fasta_panel(fasta)
        >>> p.ref_row
        ['A', 'C', 'G', 'T']
        >>> p.total_cols
        4
        >>> len(p.seq_rows)
        2
        >>> p.seq_rows[0]
        ('read1', ['A', 'C', 'T', 'T'], False)
        >>> p2 = fasta_panel(fasta, col_start=2, col_end=3)
        >>> p2.ref_row
        ['C', 'G']
    """
    seqs = read_fasta(path)
    if not seqs:
        raise ValueError(f"No sequences in {path}")

    _ref_name, ref_seq = seqs[0]

    # Slice columns if requested (1-based inclusive)
    if col_start is not None or col_end is not None:
        cs = (col_start or 1) - 1
        ce = col_end or len(ref_seq)
        ref_seq = ref_seq[cs:ce]
        seqs = [(n, s[cs:ce]) for n, s in seqs]

    aln_len = len(ref_seq)
    ref_row = list(ref_seq.upper())

    seq_rows: list[tuple[str, list[str], bool]] = []
    for name, seq in seqs[1:]:
        row = list(seq.upper()[:aln_len])
        row += ["-"] * (aln_len - len(row))
        seq_rows.append((name, row, False))

    # Column labels: 1-based position in the reference (skip gap columns)
    col_labels: list[tuple[int, str]] = []
    ref_pos = 0
    for i, base in enumerate(ref_row):
        if base != "-":
            ref_pos += 1
            if ref_pos == 1 or ref_pos % 10 == 0:
                col_labels.append((i, str(ref_pos)))

    label = Path(path).stem
    return Panel(label, ref_row, seq_rows, aln_len, col_labels)
