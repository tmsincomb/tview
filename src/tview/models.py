"""Data structures for tview alignment panels."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Panel:
    """One horizontal alignment block: a reference row + read/sequence rows.

    Attributes:
        label: Display name for the panel (e.g. filename stem).
        ref_row: Reference sequence as a list of single-character strings.
        seq_rows: Read/sequence rows as (name, bases, is_reverse) tuples.
        total_cols: Total number of display columns including insertion columns.
        col_labels: Tick positions and labels for the x-axis as (column_index, label) pairs.
        ins_columns: Column indices that represent insertion positions.

    Examples:
        >>> p = Panel("test", ["A", "C"], [("r1", ["A", "T"], False)], 2, [(0, "1")])
        >>> p.label
        'test'
        >>> p.ins_columns
        set()
    """

    label: str
    ref_row: list[str]
    seq_rows: list[tuple[str, list[str], bool]]
    total_cols: int
    col_labels: list[tuple[int, str]]
    ins_columns: set[int] | None = None

    def __post_init__(self) -> None:
        if self.ins_columns is None:
            self.ins_columns = set()
