"""Data structures for tview alignment panels.

Examples:
    >>> p = Panel("demo", ["A", "C", "G"], [], 3, [(0, "1")])
    >>> p.label
    'demo'
    >>> p.ins_columns
    set()
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Panel:
    """One horizontal alignment block: a reference row + read/sequence rows.

    Attributes:
        label: Display name for the panel (e.g. filename stem).
        ref_row: Variant-call reference sequence as single-character strings.
            Mismatches in ``seq_rows`` are computed against this row.
        seq_rows: Read/sequence rows as (name, bases, is_reverse) tuples.
        total_cols: Total number of display columns including insertion columns.
        col_labels: Tick positions and labels for the variant-call reference's
            x-axis as (column_index, label) pairs.
        ins_columns: Column indices that represent insertion positions.
        secondary_ref_row: Optional reference for orange "heterologous" highlighting.
            When a base in ``seq_rows`` mismatches ``ref_row`` but matches this row,
            it is rendered in ``heterologous_color``.
        heterologous_color: Color used for heterologous matches.
        numbering_ref_row: Optional second reference whose non-gap positions provide
            x-axis labels (e.g., HxB2 for HIV). When set, this row is rendered as
            an extra row above ``ref_row`` and its own x-axis is drawn.
        numbering_ref_label: Display name for ``numbering_ref_row`` (e.g., "HxB2").
        numbering_col_labels: Tick positions and labels for the numbering reference.
        variant_ref_label: Display name for ``ref_row`` when in dual-ref mode
            (e.g., "SF162p3_ref"). Used by ``show_row_labels`` and right-margin
            axis annotations.

    Examples:
        >>> p = Panel("test", ["A", "C"], [("r1", ["A", "T"], False)], 2, [(0, "1")])
        >>> p.label
        'test'
        >>> p.ins_columns
        set()
        >>> p.numbering_ref_row is None
        True
    """

    label: str
    ref_row: list[str]
    seq_rows: list[tuple[str, list[str], bool]]
    total_cols: int
    col_labels: list[tuple[int, str]]
    ins_columns: set[int] = field(default_factory=set)
    secondary_ref_row: list[str] | None = None
    heterologous_color: str = "#FF6F00"
    numbering_ref_row: list[str] | None = None
    numbering_ref_label: str | None = None
    numbering_col_labels: list[tuple[int, str]] | None = None
    variant_ref_label: str | None = None
