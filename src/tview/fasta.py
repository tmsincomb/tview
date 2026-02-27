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
    columns: list[int] | None = None,
) -> Panel:
    """Build a Panel from an aligned FASTA where the first sequence is the reference.

    Args:
        path: Path to the aligned FASTA file.
        columns: Sorted list of 1-based alignment column positions to include.
            Supports discrete positions, contiguous ranges, or any mix.
            When ``None``, all columns are included.

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
        >>> p2 = fasta_panel(fasta, columns=[2, 3])
        >>> p2.ref_row
        ['C', 'G']
    """
    seqs = read_fasta(path)
    if not seqs:
        raise ValueError(f"No sequences in {path}")

    _ref_name, ref_seq = seqs[0]

    # Select columns if requested (1-based positions → 0-based indices)
    orig_positions: list[int] | None = None
    if columns is not None:
        indices = sorted(i - 1 for i in columns if 1 <= i <= len(ref_seq))
        ref_seq = "".join(ref_seq[i] for i in indices)
        seqs = [(n, "".join(s[i] for i in indices if i < len(s))) for n, s in seqs]
        orig_positions = [i + 1 for i in indices]

    aln_len = len(ref_seq)
    ref_row = list(ref_seq.upper())

    seq_rows: list[tuple[str, list[str], bool]] = []
    for name, seq in seqs[1:]:
        row = list(seq.upper()[:aln_len])
        row += ["-"] * (aln_len - len(row))
        seq_rows.append((name, row, False))

    # Column labels
    col_labels: list[tuple[int, str]] = []
    if orig_positions is not None:
        # Label with original alignment positions
        for i, pos in enumerate(orig_positions):
            if i == 0 or pos % 10 == 0:
                col_labels.append((i, str(pos)))
    else:
        # Default: 1-based reference position (skip gap columns)
        ref_pos = 0
        for i, base in enumerate(ref_row):
            if base != "-":
                ref_pos += 1
                if ref_pos == 1 or ref_pos % 10 == 0:
                    col_labels.append((i, str(ref_pos)))

    label = Path(path).stem
    return Panel(label, ref_row, seq_rows, aln_len, col_labels)
