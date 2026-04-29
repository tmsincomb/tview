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


def _build_ref_col_labels(
    ref_row: list[str], tick_every: int = 10
) -> list[tuple[int, str]]:
    """Walk a reference row's non-gap positions, labeling on a tick interval.

    Args:
        ref_row: Reference sequence as single-character strings.
        tick_every: Label every Nth non-gap position (plus position 1).
            ``1`` labels every column; ``10`` labels every tenth.

    Returns:
        List of (column_index, position_label) tuples.
    """
    labels: list[tuple[int, str]] = []
    ref_pos = 0
    for i, base in enumerate(ref_row):
        if base != "-":
            ref_pos += 1
            if ref_pos == 1 or ref_pos % tick_every == 0:
                labels.append((i, str(ref_pos)))
    return labels


def fasta_panel(
    path: str | Path,
    columns: list[int] | None = None,
    max_rows: int | None = None,
    variant_ref: str | None = None,
    numbering_ref: str | None = None,
    tick_every: int = 10,
) -> Panel:
    """Build a Panel from an aligned FASTA.

    By default the first sequence is treated as the variant-call reference and
    also drives the x-axis numbering. Pass ``variant_ref`` to pick a different
    sequence by FASTA header name, and ``numbering_ref`` to layer a second
    reference whose non-gap positions provide an additional x-axis (typical
    HIV use case: HxB2 numbering with variant calls against a different strain).

    Args:
        path: Path to the aligned FASTA file.
        columns: Sorted list of 1-based alignment column positions to include.
            Supports discrete positions, contiguous ranges, or any mix.
            When ``None``, all columns are included.
        max_rows: Maximum number of non-reference rows to keep. References
            (variant + numbering) are always retained. When dual-ref mode is
            active, the cap counts samples only. In single-ref mode (no
            ``variant_ref``/``numbering_ref`` given), the cap counts the
            reference plus samples for backward compatibility.
        variant_ref: FASTA header name for the variant-call reference. Mismatches
            are computed against this row. Defaults to the first sequence.
        numbering_ref: FASTA header name for the numbering reference. When set
            and different from ``variant_ref``, this row is rendered as an extra
            row above ``ref_row`` and an additional x-axis is drawn.

    Returns:
        A Panel with reference row, sequence rows, and column labels. When dual
        refs are active, ``numbering_ref_row`` and ``numbering_col_labels`` are
        also populated.

    Raises:
        ValueError: If the FASTA contains no sequences, or a referenced header
            name is not found.

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
        >>> p.numbering_ref_row is None
        True
        >>> p2 = fasta_panel(fasta, columns=[2, 3])
        >>> p2.ref_row
        ['C', 'G']
    """
    seqs = read_fasta(path)
    if not seqs:
        raise ValueError(f"No sequences in {path}")

    names = [name for name, _ in seqs]

    # Resolve variant-ref index (defaults to 0).
    if variant_ref is None:
        variant_ref_idx = 0
    else:
        try:
            variant_ref_idx = names.index(variant_ref)
        except ValueError:
            raise ValueError(
                f"variant-ref {variant_ref!r} not found in {path}. "
                f"Available headers: {', '.join(names)}"
            ) from None

    # Resolve numbering-ref index. None or matching the variant ref → single-ref mode.
    numbering_ref_idx: int | None = None
    if numbering_ref is not None and numbering_ref != names[variant_ref_idx]:
        try:
            numbering_ref_idx = names.index(numbering_ref)
        except ValueError:
            raise ValueError(
                f"numbering-ref {numbering_ref!r} not found in {path}. "
                f"Available headers: {', '.join(names)}"
            ) from None
        if numbering_ref_idx == variant_ref_idx:
            numbering_ref_idx = None

    dual_ref = numbering_ref_idx is not None

    # Truncate sample pool by max_rows.
    if dual_ref or variant_ref is not None:
        # New code path: refs are not counted toward max_rows.
        excluded = {variant_ref_idx}
        if numbering_ref_idx is not None:
            excluded.add(numbering_ref_idx)
        pool: list[tuple[str, str]] = [
            (name, seq) for i, (name, seq) in enumerate(seqs) if i not in excluded
        ]
        if max_rows is not None and max_rows >= 0:
            pool = pool[:max_rows]
    else:
        # Backward-compat path: preserve existing seqs[: max_rows + 1] semantics
        # exactly when no new flags are supplied.
        if max_rows is not None and max_rows >= 0:
            seqs = seqs[: max_rows + 1]
        pool = [(name, seq) for i, (name, seq) in enumerate(seqs) if i != 0]

    variant_ref_name = names[variant_ref_idx]
    variant_seq = seqs[variant_ref_idx][1]
    numbering_ref_name: str | None = None
    numbering_seq: str | None = None
    if numbering_ref_idx is not None:
        numbering_ref_name = names[numbering_ref_idx]
        numbering_seq = seqs[numbering_ref_idx][1]

    # Select columns if requested (1-based positions → 0-based indices).
    orig_positions: list[int] | None = None
    if columns is not None:
        indices = sorted(i - 1 for i in columns if 1 <= i <= len(variant_seq))
        variant_seq = "".join(variant_seq[i] for i in indices)
        if numbering_seq is not None:
            numbering_seq = "".join(
                numbering_seq[i] for i in indices if i < len(numbering_seq)
            )
        pool = [(n, "".join(s[i] for i in indices if i < len(s))) for n, s in pool]
        orig_positions = [i + 1 for i in indices]

    aln_len = len(variant_seq)
    ref_row = list(variant_seq.upper())

    seq_rows: list[tuple[str, list[str], bool]] = []
    for name, seq in pool:
        row = list(seq.upper()[:aln_len])
        row += ["-"] * (aln_len - len(row))
        seq_rows.append((name, row, False))

    # Variant-ref column labels.
    col_labels: list[tuple[int, str]] = []
    if orig_positions is not None:
        for i, pos in enumerate(orig_positions):
            if i == 0 or pos % tick_every == 0:
                col_labels.append((i, str(pos)))
    else:
        col_labels = _build_ref_col_labels(ref_row, tick_every=tick_every)

    # Numbering-ref row + labels (always by non-gap positions in numbering ref).
    numbering_ref_row: list[str] | None = None
    numbering_col_labels: list[tuple[int, str]] | None = None
    numbering_ref_label: str | None = None
    variant_ref_label: str | None = None
    if numbering_seq is not None and numbering_ref_name is not None:
        numbering_ref_row = list(numbering_seq.upper())
        numbering_ref_row += ["-"] * (aln_len - len(numbering_ref_row))
        numbering_ref_row = numbering_ref_row[:aln_len]
        numbering_col_labels = _build_ref_col_labels(
            numbering_ref_row, tick_every=tick_every
        )
        numbering_ref_label = numbering_ref_name
        variant_ref_label = variant_ref_name

    label = Path(path).stem
    return Panel(
        label=label,
        ref_row=ref_row,
        seq_rows=seq_rows,
        total_cols=aln_len,
        col_labels=col_labels,
        numbering_ref_row=numbering_ref_row,
        numbering_ref_label=numbering_ref_label,
        numbering_col_labels=numbering_col_labels,
        variant_ref_label=variant_ref_label,
    )
