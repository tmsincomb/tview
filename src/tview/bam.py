"""BAM parsing and panel construction.

Provides CIGAR-aware alignment extraction and panel building from indexed
BAM files.  Requires ``pysam`` at runtime (imported lazily in ``bam_panel``).

Examples:
    >>> CIGAR_MATCH, CIGAR_INS, CIGAR_DEL
    (0, 1, 2)
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

from tview.models import Panel

# -- CIGAR operations --------------------------------------------------
CIGAR_MATCH = 0  # M
CIGAR_INS = 1  # I
CIGAR_DEL = 2  # D
CIGAR_REF_SKIP = 3  # N
CIGAR_SOFT_CLIP = 4  # S
CIGAR_SEQ_MATCH = 7  # =
CIGAR_SEQ_MISMATCH = 8  # X


def build_read_row(
    read: Any,
    ref_start: int,
    ref_end: int,
) -> tuple[dict[int, str], dict[int, list[str]]]:
    """Extract aligned bases and insertions from a single pysam read.

    Walks the CIGAR string to map query bases onto reference positions,
    collecting insertions keyed by their anchor reference position.

    Args:
        read: A pysam.AlignedSegment with cigartuples and query_sequence.
        ref_start: 0-based start of the reference window (inclusive).
        ref_end: 0-based end of the reference window (exclusive).

    Returns:
        A tuple of (aligned, inserts) where aligned maps ref positions to
        bases and inserts maps ref positions to lists of inserted bases.
    """
    aligned: dict[int, str] = {}
    inserts: dict[int, list[str]] = defaultdict(list)
    qpos, rpos = 0, read.reference_start
    for op, length in read.cigartuples:
        if op in (CIGAR_MATCH, CIGAR_SEQ_MATCH, CIGAR_SEQ_MISMATCH):
            for _ in range(length):
                if ref_start <= rpos < ref_end:
                    aligned[rpos] = read.query_sequence[qpos].upper()
                qpos += 1
                rpos += 1
        elif op == CIGAR_INS:
            anchor = rpos - 1
            if ref_start <= anchor < ref_end:
                for j in range(length):
                    inserts[anchor].append(read.query_sequence[qpos + j].upper())
            qpos += length
        elif op == CIGAR_DEL:
            for _ in range(length):
                if ref_start <= rpos < ref_end:
                    aligned[rpos] = "-"
                rpos += 1
        elif op == CIGAR_REF_SKIP:
            rpos += length
        elif op == CIGAR_SOFT_CLIP:
            qpos += length
    return aligned, inserts


def bam_panel(bam_path: str | Path, ref_path: str | Path, region: str) -> Panel:
    """Build a Panel from a BAM file with reference FASTA and genomic region.

    Reads are sorted by start position and strand. Insertion columns are
    expanded so all reads align on a common grid.

    Args:
        bam_path: Path to the indexed BAM file.
        ref_path: Path to the reference FASTA (must be indexed).
        region: Genomic region string in "chrom:start-end" format (0-based start).

    Returns:
        A Panel with reference row, read rows, insertion columns, and tick labels.
    """
    import pysam

    chrom, rest = region.split(":")
    start, end = [int(x) for x in rest.replace(",", "").split("-")]

    with pysam.FastaFile(ref_path) as fasta:
        ref_seq = fasta.fetch(chrom, start, end).upper()

    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        reads = [
            r
            for r in samfile.fetch(chrom, start, end)
            if not r.is_unmapped and r.cigartuples
        ]
    reads.sort(key=lambda r: (r.reference_start, r.is_reverse))

    # Find max insertion at each ref position
    max_ins: dict[int, int] = defaultdict(int)
    read_data = []
    for read in reads:
        aligned, inserts = build_read_row(read, start, end)
        read_data.append((read, aligned, inserts))
        for rpos, bases in inserts.items():
            max_ins[rpos] = max(max_ins[rpos], len(bases))

    # Build column map
    col_map: dict[int, int] = {}
    ins_col_set: set[int] = set()
    col = 0
    for rpos in range(start, end):
        col_map[rpos] = col
        col += 1
        n_ins = max_ins.get(rpos, 0)
        for j in range(n_ins):
            ins_col_set.add(col + j)
        col += n_ins
    total_cols = col

    # Build ref row with '-' in insertion columns
    ref_row: list[str] = []
    for rpos in range(start, end):
        ref_row.append(ref_seq[rpos - start])
        for _ in range(max_ins.get(rpos, 0)):
            ref_row.append("-")

    # Build sequence rows
    seq_rows: list[tuple[str, list[str], bool]] = []
    for read, aligned, inserts in read_data:
        row = [" "] * total_cols
        for rpos in range(start, end):
            c = col_map[rpos]
            if rpos in aligned:
                row[c] = aligned[rpos]
            if rpos in aligned or rpos in inserts:
                n_ins = max_ins.get(rpos, 0)
                read_ins = inserts.get(rpos, [])
                for j in range(n_ins):
                    if j < len(read_ins):
                        row[c + 1 + j] = read_ins[j]
                    else:
                        row[c + 1 + j] = "-"
        seq_rows.append((read.query_name, row, read.is_reverse))

    # Column labels: 1-based relative, ticks at 1, 10, 20...
    ref_width = end - start
    tick_1based = [1] + list(range(10, ref_width + 1, 10))
    col_labels = [
        (col_map[start + p - 1], str(p)) for p in tick_1based if (start + p - 1) < end
    ]

    label = Path(bam_path).stem
    return Panel(label, ref_row, seq_rows, total_cols, col_labels, ins_col_set)
