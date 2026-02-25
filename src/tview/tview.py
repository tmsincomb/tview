#!/usr/bin/env python3
"""
tview_image.py -- Publication-quality alignment viewer.

Supports BAM files (with ref FASTA) and pre-aligned FASTA (e.g. MAFFT output).
Multiple inputs can be stacked vertically in a single figure.

Usage:
  # Single BAM
  python tview_image.py --bam sample.bam --ref ref.fa --region chr1:1-50 -o out.png

  # Stacked BAMs (shared ref + region)
  python tview_image.py --bam a.bam b.bam --ref ref.fa --region chr1:1-50 -o out.png

  # Aligned FASTA (first sequence = reference)
  python tview_image.py --fasta aligned.fasta -o out.png --palette aa

  # Stacked FASTAs
  python tview_image.py --fasta group1.fasta group2.fasta -o out.png --palette aa

  # FASTA with column range (1-based, inclusive)
  python tview_image.py --fasta aligned.fasta --columns 1-120 -o out.png

  # Mix BAM and FASTA panels
  python tview_image.py --bam a.bam --ref ref.fa --region chr1:1-50 --fasta aln.fasta -o out.png
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

# -- Color schemes -----------------------------------------------------
NT_COLORS = {
    "A": "#4CAF50",
    "C": "#2196F3",
    "G": "#FF9800",
    "T": "#F44336",
    "N": "#9E9E9E",
    "-": "#9E9E9E",
}
AA_COLORS = {
    "A": "#2196F3",
    "V": "#2196F3",
    "L": "#2196F3",
    "I": "#2196F3",
    "M": "#2196F3",
    "F": "#2196F3",
    "W": "#2196F3",
    "P": "#2196F3",
    "K": "#F44336",
    "R": "#F44336",
    "H": "#F44336",
    "D": "#E040FB",
    "E": "#E040FB",
    "S": "#4CAF50",
    "T": "#4CAF50",
    "N": "#4CAF50",
    "Q": "#4CAF50",
    "G": "#FF9800",
    "C": "#FF9800",
    "Y": "#FF9800",
    "*": "#9E9E9E",
    "-": "#9E9E9E",
    "X": "#9E9E9E",
}
MISMATCH_BG = "#FFEB3B55"
INS_BG = "#CE93D833"
GAP_COLOR = "#9E9E9E"
FWD_ALPHA = 1.0
REV_ALPHA = 0.85
SEPARATOR_COLOR = "#BDBDBD"

# -- CIGAR operations --------------------------------------------------
CIGAR_MATCH = 0  # M
CIGAR_INS = 1  # I
CIGAR_DEL = 2  # D
CIGAR_REF_SKIP = 3  # N
CIGAR_SOFT_CLIP = 4  # S
CIGAR_SEQ_MATCH = 7  # =
CIGAR_SEQ_MISMATCH = 8  # X


# ======================================================================
#  Data structures -- each input becomes a Panel
# ======================================================================
@dataclass
class Panel:
    """One horizontal alignment block: a reference row + read/sequence rows."""

    label: str
    ref_row: list[str]
    seq_rows: list[tuple[str, list[str], bool]]
    total_cols: int
    col_labels: list[tuple[int, str]]
    ins_columns: set[int] | None = None

    def __post_init__(self) -> None:
        if self.ins_columns is None:
            self.ins_columns = set()


# ======================================================================
#  FASTA panel builder
# ======================================================================
def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Parse FASTA -> list[(name, sequence)]."""
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
    """Build a Panel from an aligned FASTA. First sequence = reference.

    col_start / col_end are 1-based inclusive column indices into the
    alignment (after MAFFT gaps).
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


# ======================================================================
#  BAM panel builder
# ======================================================================
def build_read_row(
    read,
    ref_start: int,
    ref_end: int,
) -> tuple[dict[int, str], dict[int, list[str]]]:
    """Extract aligned bases and insertions from a single read."""
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
    """Build a Panel from a BAM file with reference FASTA and genomic region."""
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


# ======================================================================
#  Renderer
# ======================================================================
def _resolve_font(
    fontsize: float,
) -> tuple[fm.FontProperties, fm.FontProperties]:
    """Resolve font, preferring Helvetica (macOS), falling back to DejaVu Sans Mono."""
    helv_path = fm.findfont(fm.FontProperties(family="Helvetica", style="normal"))
    if "Helvetica" in helv_path:
        mono = fm.FontProperties(fname=helv_path, size=fontsize, weight="bold")
        mono_sm = fm.FontProperties(fname=helv_path, size=fontsize, weight="bold")
        return mono, mono_sm

    mono_bold_path = fm.findfont(
        fm.FontProperties(family="monospace", weight="bold", style="normal")
    )
    if "Oblique" in mono_bold_path or "Italic" in mono_bold_path:
        import glob as _gl

        candidates = _gl.glob(
            str(Path(fm.findfont("monospace")).parent / "DejaVuSansMono-Bold.ttf")
        )
        if candidates:
            mono_bold_path = candidates[0]

    mono = fm.FontProperties(fname=mono_bold_path, size=fontsize)
    mono_sm = fm.FontProperties(fname=mono_bold_path, size=fontsize)
    return mono, mono_sm


def render_panels(
    panels: list[Panel],
    out_path: str | Path = "alignment.png",
    fontsize: float = 12,
    dpi: int = 600,
    palette: str = "nt",
    cell: float | None = None,
) -> None:
    """Render alignment panels to an image file."""
    if cell is None:
        cell = fontsize / 72  # 1 pt = 1/72 inch -> cell fits one character
    colors = AA_COLORS if palette == "aa" else NT_COLORS
    mono, mono_sm = _resolve_font(fontsize)

    # Compute total height: for each panel, 1 ref row + N seq rows + separator
    max_cols = max(p.total_cols for p in panels)
    total_rows = 0
    panel_y_offsets: list[int] = []
    for i, p in enumerate(panels):
        panel_y_offsets.append(total_rows)
        total_rows += 1 + len(p.seq_rows)
        if i < len(panels) - 1:
            total_rows += 1

    fig_w = max(4, max_cols * cell + 0.5)
    fig_h = max(1.0, total_rows * cell + 0.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-0.5, max_cols - 0.5)
    ax.set_ylim(total_rows - 0.5, -0.5)
    ax.set_aspect("equal")

    for pi, panel in enumerate(panels):
        y0 = panel_y_offsets[pi]
        n_panel_rows = 1 + len(panel.seq_rows)

        # Shade insertion columns
        for ic in panel.ins_columns:
            ax.add_patch(
                plt.Rectangle(
                    (ic - 0.5, y0 - 0.5),
                    1,
                    n_panel_rows,
                    facecolor=INS_BG,
                    edgecolor="none",
                    zorder=0,
                )
            )

        # Reference row
        for c, base in enumerate(panel.ref_row):
            clr = "#000000" if base == "-" else colors.get(base, "#9E9E9E")
            ax.text(
                c, y0, base, ha="center", va="center", fontproperties=mono, color=clr
            )

        # Sequence rows
        for ri, (name, row, is_rev) in enumerate(panel.seq_rows):
            y = y0 + 1 + ri
            alpha = REV_ALPHA if is_rev else FWD_ALPHA
            strand_char = "," if is_rev else "."

            for c, base in enumerate(row):
                if base == " ":
                    continue
                ref_base = panel.ref_row[c] if c < len(panel.ref_row) else "-"

                if base == "-":
                    ax.text(
                        c,
                        y,
                        "-",
                        ha="center",
                        va="center",
                        fontproperties=mono,
                        color="#000000",
                        alpha=alpha,
                    )
                elif base == ref_base:
                    ax.text(
                        c,
                        y,
                        strand_char,
                        ha="center",
                        va="center",
                        fontproperties=mono,
                        color="#000000",
                        alpha=alpha,
                    )
                else:
                    ax.add_patch(
                        plt.Rectangle(
                            (c - 0.5, y - 0.5),
                            1,
                            1,
                            facecolor=MISMATCH_BG,
                            edgecolor="none",
                        )
                    )
                    display = base.lower() if is_rev else base
                    ax.text(
                        c,
                        y,
                        display,
                        ha="center",
                        va="center",
                        fontproperties=mono,
                        color=colors.get(base, "#000000"),
                        alpha=alpha,
                    )

        # Panel label (left side)
        if len(panels) > 1:
            ax.text(
                -1.5,
                y0 + n_panel_rows / 2 - 0.5,
                panel.label,
                ha="right",
                va="center",
                fontproperties=mono,
                color="#616161",
            )

        # Separator line between panels
        if pi < len(panels) - 1:
            sep_y = y0 + n_panel_rows + 0.0
            ax.axhline(y=sep_y, color=SEPARATOR_COLOR, lw=0.5, ls="-", xmin=0, xmax=1)

    # X-axis from first panel, placed on top
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    first = panels[0]
    tick_idx = [ci for ci, _ in first.col_labels]
    tick_lbl = [lb for _, lb in first.col_labels]
    ax.set_xticks(tick_idx)
    ax.set_xticklabels(tick_lbl, rotation=0, ha="center", fontproperties=mono_sm)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.01)
    plt.savefig(
        out_path,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0.05,
        facecolor="white",
        transparent=False,
    )
    plt.close()
    print(f"Saved: {out_path} ({dpi} dpi, {len(panels)} panel(s), " f"{max_cols} cols)")


if __name__ == "__main__":
    from tview.cli import main

    main()
