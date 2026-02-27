"""Rendering engine for alignment panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from tview.config import (
    AA_COLORS,
    FALLBACK_BASE_COLOR,
    FONT_FALLBACK_FILENAME,
    FONT_PREFERENCES,
    FWD_ALPHA,
    INS_BG,
    MISMATCH_BG,
    NT_COLORS,
    PANEL_LABEL_COLOR,
    REV_ALPHA,
    SEPARATOR_COLOR,
    TEXT_COLOR,
)
from tview.models import Panel


def _resolve_font(
    fontsize: float,
) -> tuple[fm.FontProperties, fm.FontProperties]:
    """Resolve monospace font for alignment rendering.

    Tries each font in ``FONT_PREFERENCES`` (from style.yaml) in order, then
    falls back to ``FONT_FALLBACK_FILENAME``.

    Args:
        fontsize: Font size in points.

    Returns:
        A tuple of (mono, mono_sm) FontProperties for base text and tick labels.
    """
    for pref in FONT_PREFERENCES:
        family = pref["family"]
        weight = pref.get("weight", "normal")
        found_path = fm.findfont(fm.FontProperties(family=family, style="normal"))
        if family in found_path:
            mono = fm.FontProperties(fname=found_path, size=fontsize, weight=weight)
            mono_sm = fm.FontProperties(fname=found_path, size=fontsize, weight=weight)
            return mono, mono_sm

    # Final fallback: probe for the fallback font file
    mono_bold_path = fm.findfont(
        fm.FontProperties(family="monospace", weight="bold", style="normal")
    )
    if "Oblique" in mono_bold_path or "Italic" in mono_bold_path:
        fallback = Path(fm.findfont("monospace")).parent / FONT_FALLBACK_FILENAME
        if fallback.exists():
            mono_bold_path = str(fallback)

    mono = fm.FontProperties(fname=mono_bold_path, size=fontsize)
    mono_sm = fm.FontProperties(fname=mono_bold_path, size=fontsize)
    return mono, mono_sm


def panel_figsize(
    panels: list[Panel],
    fontsize: float = 12,
    cell: float | None = None,
) -> tuple[float, float]:
    """Compute recommended figure size for a set of alignment panels.

    Useful when creating external figures or patchworklib Bricks that need
    to match the natural size of the alignment rendering.

    Args:
        panels: List of Panel objects to measure.
        fontsize: Font size in points for base characters.
        cell: Cell size in inches. Defaults to fontsize / 72.

    Returns:
        A (width, height) tuple in inches.

    Examples:
        >>> from tview.models import Panel
        >>> p = Panel("t", list("ACGT"), [("s1", list("ACGT"), False)], 4, [(0, "1")])
        >>> w, h = panel_figsize([p])
        >>> w >= 4
        True
    """
    if cell is None:
        cell = fontsize / 72
    max_cols = max(p.total_cols for p in panels)
    total_rows = 0
    for i, p in enumerate(panels):
        total_rows += 1 + len(p.seq_rows)
        if i < len(panels) - 1:
            total_rows += 1
    fig_w = max(4, max_cols * cell + 0.5)
    fig_h = max(1.0, total_rows * cell + 0.6)
    return (fig_w, fig_h)


def draw_panels(
    panels: list[Panel],
    ax: Axes,
    fontsize: float = 12,
    palette: str = "nt",
    cell: float | None = None,
    classic: bool = False,
) -> Axes:
    """Draw alignment panels onto the given axes.

    Draws reference rows, sequence rows, mismatch highlights, insertion
    column shading, panel labels, separator lines, and tick configuration
    onto *ax*. The caller is responsible for figure creation and saving.

    Compatible with any ``matplotlib.axes.Axes`` subclass, including
    ``patchworklib.Brick`` objects and standard subplot axes.

    Args:
        panels: List of Panel objects to render vertically.
        ax: Matplotlib axes (or compatible subclass) to draw on.
        fontsize: Font size in points for base characters.
        palette: Color scheme, either ``"nt"`` for nucleotides or ``"aa"`` for amino acids.
        cell: Cell size in inches. Defaults to fontsize / 72.
        classic: When True, render in black-and-white with no color highlighting.

    Returns:
        The axes object (same as *ax*), for method chaining.

    Examples:
        >>> import matplotlib
        >>> matplotlib.use("Agg")
        >>> import matplotlib.pyplot as plt
        >>> from tview.models import Panel
        >>> fig, ax = plt.subplots(figsize=(4, 1))
        >>> p = Panel("t", list("ACGT"), [("s1", list("ACGT"), False)], 4, [(0, "1")])
        >>> result = draw_panels([p], ax)
        >>> result is ax
        True
        >>> plt.close(fig)
    """
    if cell is None:
        cell = fontsize / 72
    colors = AA_COLORS if palette == "aa" else NT_COLORS

    if classic:
        colors = {k: "#000000" for k in colors}
        mismatch_bg = "#FFFFFF00"
        ins_bg = "#FFFFFF00"
    else:
        mismatch_bg = MISMATCH_BG
        ins_bg = INS_BG
    mono, mono_sm = _resolve_font(fontsize)

    max_cols = max(p.total_cols for p in panels)
    total_rows = 0
    panel_y_offsets: list[int] = []
    for i, p in enumerate(panels):
        panel_y_offsets.append(total_rows)
        total_rows += 1 + len(p.seq_rows)
        if i < len(panels) - 1:
            total_rows += 1

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
                    facecolor=ins_bg,
                    edgecolor="none",
                    zorder=0,
                )
            )

        # Reference row
        for c, base in enumerate(panel.ref_row):
            clr = TEXT_COLOR if base == "-" else colors.get(base, FALLBACK_BASE_COLOR)
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
                        color=TEXT_COLOR,
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
                        color=TEXT_COLOR,
                        alpha=alpha,
                    )
                else:
                    ax.add_patch(
                        plt.Rectangle(
                            (c - 0.5, y - 0.5),
                            1,
                            1,
                            facecolor=mismatch_bg,
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
                        color=colors.get(base, TEXT_COLOR),
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
                color=PANEL_LABEL_COLOR,
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

    return ax


def render_panels(
    panels: list[Panel],
    out_path: str | Path = "alignment.png",
    fontsize: float = 12,
    dpi: int = 600,
    palette: str = "nt",
    cell: float | None = None,
    classic: bool = False,
) -> None:
    """Render alignment panels to a publication-quality image file.

    Convenience wrapper around :func:`draw_panels` that handles figure
    creation, layout adjustment, and file saving.

    Each panel is drawn as a reference row followed by read rows. Matches
    are shown as dots (forward) or commas (reverse), mismatches are
    highlighted with colored backgrounds, and insertion columns are shaded.

    Args:
        panels: List of Panel objects to render vertically.
        out_path: Output image path (format inferred from extension).
        fontsize: Font size in points for base characters.
        dpi: Output resolution in dots per inch.
        palette: Color scheme, either ``"nt"`` for nucleotides or ``"aa"`` for amino acids.
        cell: Cell size in inches. Defaults to fontsize / 72.
        classic: When True, render in black-and-white with no color highlighting.
    """
    fig_w, fig_h = panel_figsize(panels, fontsize, cell)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    draw_panels(
        panels, ax, fontsize=fontsize, palette=palette, cell=cell, classic=classic
    )
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
    max_cols = max(p.total_cols for p in panels)
    print(f"Saved: {out_path} ({dpi} dpi, {len(panels)} panel(s), {max_cols} cols)")
