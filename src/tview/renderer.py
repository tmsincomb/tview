"""Rendering engine for alignment panels."""

from __future__ import annotations

import logging
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
    GAP_COLOR,
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
            mono_sm = fm.FontProperties(
                fname=found_path, size=fontsize * 0.8, weight=weight
            )
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
    mono_sm = fm.FontProperties(fname=mono_bold_path, size=fontsize * 0.8)
    return mono, mono_sm


def _per_panel_axis_rows(panel: Panel, per_panel_axes: bool) -> tuple[int, int]:
    """Return (top_axis_rows, bottom_axis_rows) added per panel when per-panel
    axes are enabled.

    A bottom inline tick row is only added for dual-ref panels (where
    ``numbering_col_labels`` is populated). Single-ref panels get a top tick
    row only, matching the existing single-axis-on-top default.
    """
    if not per_panel_axes:
        return (0, 0)
    has_bottom = panel.numbering_col_labels is not None
    return (1, 1 if has_bottom else 0)


def panel_figsize(
    panels: list[Panel],
    fontsize: float = 12,
    cell: float | None = None,
    per_panel_axes: bool = False,
) -> tuple[float, float]:
    """Compute recommended figure size for a set of alignment panels.

    Useful when creating external figures or patchworklib Bricks that need
    to match the natural size of the alignment rendering.

    Args:
        panels: List of Panel objects to measure.
        fontsize: Font size in points for base characters.
        cell: Cell size in inches. Defaults to fontsize / 72.
        per_panel_axes: When True, reserve extra rows for inline tick labels
            above (and, for dual-ref panels, below) each panel.

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
        n_top_refs = 1 + (1 if p.numbering_ref_row is not None else 0)
        n_axis_top, n_axis_bot = _per_panel_axis_rows(p, per_panel_axes)
        total_rows += n_axis_top + n_top_refs + len(p.seq_rows) + n_axis_bot
        if i < len(panels) - 1:
            total_rows += 1
    has_numbering = any(p.numbering_ref_row is not None for p in panels)
    extra_axis_pad = 0.0 if per_panel_axes else (0.4 if has_numbering else 0.0)
    fig_w = max(4, max_cols * cell + 0.5)
    fig_h = max(1.0, total_rows * cell + 0.6 + extra_axis_pad)
    return (fig_w, fig_h)


def draw_panels(
    panels: list[Panel],
    ax: Axes,
    fontsize: float = 12,
    palette: str = "nt",
    cell: float | None = None,
    classic: bool = False,
    show_row_labels: bool = False,
    row_label_width: float = 8.0,
    per_panel_axes: bool = False,
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
        per_panel_axes: When True, draw an inline tick-label row above each
            panel's top reference (and, for dual-ref panels, below each panel's
            last sample row) instead of placing tick labels only at the figure
            edges.

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
    panel_n_top_refs: list[int] = []
    panel_n_axis_top: list[int] = []
    panel_n_axis_bot: list[int] = []
    for i, p in enumerate(panels):
        panel_y_offsets.append(total_rows)
        n_top_refs = 1 + (1 if p.numbering_ref_row is not None else 0)
        n_axis_top, n_axis_bot = _per_panel_axis_rows(p, per_panel_axes)
        panel_n_top_refs.append(n_top_refs)
        panel_n_axis_top.append(n_axis_top)
        panel_n_axis_bot.append(n_axis_bot)
        total_rows += n_axis_top + n_top_refs + len(p.seq_rows) + n_axis_bot
        if i < len(panels) - 1:
            total_rows += 1

    left_pad = -row_label_width if show_row_labels else -0.5
    ax.set_xlim(left_pad, max_cols - 0.5)
    ax.set_ylim(total_rows - 0.5, -0.5)
    ax.set_aspect("equal")

    for pi, panel in enumerate(panels):
        block_top_y = panel_y_offsets[pi]
        n_top_refs = panel_n_top_refs[pi]
        n_axis_top = panel_n_axis_top[pi]
        n_axis_bot = panel_n_axis_bot[pi]
        axis_top_y = block_top_y  # inline tick row (only used if per_panel_axes)
        y0 = block_top_y + n_axis_top  # first content row (numbering_ref or single ref)
        n_panel_rows = n_top_refs + len(panel.seq_rows)
        # y position of the variant-call ref_row (used for mismatch baseline).
        ref_y = y0 + (n_top_refs - 1)
        # y position of the inline bottom tick row (only used if dual-ref + flag on).
        axis_bot_y = ref_y + len(panel.seq_rows) + 1

        # Inline per-panel tick labels (only when per_panel_axes is on).
        if per_panel_axes:
            top_labels = panel.numbering_col_labels or panel.col_labels
            for ci, lbl in top_labels:
                ax.text(
                    ci,
                    axis_top_y,
                    lbl,
                    ha="center",
                    va="center",
                    fontproperties=mono_sm,
                    color=PANEL_LABEL_COLOR,
                )
            if n_axis_bot:
                for ci, lbl in panel.col_labels:
                    ax.text(
                        ci,
                        axis_bot_y,
                        lbl,
                        ha="center",
                        va="center",
                        fontproperties=mono_sm,
                        color=PANEL_LABEL_COLOR,
                    )
            # Right-margin annotations naming the axes (one set per panel).
            if panel.numbering_ref_label:
                ax.text(
                    max_cols - 0.2,
                    axis_top_y,
                    f"({panel.numbering_ref_label})",
                    ha="left",
                    va="center",
                    fontproperties=mono_sm,
                    color=PANEL_LABEL_COLOR,
                    fontstyle="italic",
                )
            if n_axis_bot and panel.variant_ref_label:
                ax.text(
                    max_cols - 0.2,
                    axis_bot_y,
                    f"({panel.variant_ref_label})",
                    ha="left",
                    va="center",
                    fontproperties=mono_sm,
                    color=PANEL_LABEL_COLOR,
                    fontstyle="italic",
                )

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

        # Numbering reference row (no mismatch coloring — purely for numbering).
        if panel.numbering_ref_row is not None:
            for c, base in enumerate(panel.numbering_ref_row):
                clr = (
                    GAP_COLOR if base == "-" else colors.get(base, FALLBACK_BASE_COLOR)
                )
                ax.text(
                    c,
                    y0,
                    base,
                    ha="center",
                    va="center",
                    fontproperties=mono,
                    color=clr,
                )

        # Variant-call reference row.
        for c, base in enumerate(panel.ref_row):
            clr = GAP_COLOR if base == "-" else colors.get(base, FALLBACK_BASE_COLOR)
            ax.text(
                c,
                ref_y,
                base,
                ha="center",
                va="center",
                fontproperties=mono,
                color=clr,
            )

        # Sequence rows
        for ri, (name, row, is_rev) in enumerate(panel.seq_rows):
            y = ref_y + 1 + ri
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
                        color=GAP_COLOR,
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
                    sec_ref = panel.secondary_ref_row
                    is_heterologous = (
                        not classic
                        and sec_ref is not None
                        and c < len(sec_ref)
                        and base == sec_ref[c]
                    )
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
                    letter_color = (
                        panel.heterologous_color
                        if is_heterologous
                        else colors.get(base, TEXT_COLOR)
                    )
                    ax.text(
                        c,
                        y,
                        display,
                        ha="center",
                        va="center",
                        fontproperties=mono,
                        color=letter_color,
                        alpha=alpha,
                        fontweight="bold" if is_heterologous else "normal",
                    )

        # Panel label (left side) — only when not showing per-row labels
        if len(panels) > 1 and not show_row_labels:
            ax.text(
                -1.5,
                y0 + n_panel_rows / 2 - 0.5,
                panel.label,
                ha="right",
                va="center",
                fontproperties=mono,
                color=PANEL_LABEL_COLOR,
            )

        # Per-row labels on the left (one per ref + sequence row)
        if show_row_labels:
            label_x = -1.0
            if panel.numbering_ref_row is not None:
                ax.text(
                    label_x,
                    y0,
                    panel.numbering_ref_label or panel.label,
                    ha="right",
                    va="center",
                    fontproperties=mono_sm,
                    color=PANEL_LABEL_COLOR,
                    fontweight="bold",
                )
                ax.text(
                    label_x,
                    ref_y,
                    panel.variant_ref_label or panel.label,
                    ha="right",
                    va="center",
                    fontproperties=mono_sm,
                    color=PANEL_LABEL_COLOR,
                    fontweight="bold",
                )
            else:
                ax.text(
                    label_x,
                    y0,
                    panel.label,
                    ha="right",
                    va="center",
                    fontproperties=mono_sm,
                    color=PANEL_LABEL_COLOR,
                    fontweight="bold",
                )
            for ri, (name, _row, _is_rev) in enumerate(panel.seq_rows):
                ax.text(
                    label_x,
                    ref_y + 1 + ri,
                    name,
                    ha="right",
                    va="center",
                    fontproperties=mono_sm,
                    color=PANEL_LABEL_COLOR,
                )

        # Separator line between panels
        if pi < len(panels) - 1:
            sep_y = y0 + n_panel_rows + n_axis_bot
            ax.axhline(y=sep_y, color=SEPARATOR_COLOR, lw=0.5, ls="-", xmin=0, xmax=1)

    if per_panel_axes:
        # Inline per-panel ticks have already been drawn inside the loop.
        # Hide the native matplotlib axes entirely.
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(axis="x", length=0)
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax

    # X-axis configuration. In dual-ref mode, top axis = numbering ref;
    # bottom secondary axis = variant-call ref. Otherwise single axis on top.
    first = panels[0]
    has_numbering = any(p.numbering_col_labels is not None for p in panels)

    if has_numbering and first.numbering_col_labels is not None:
        top_labels = first.numbering_col_labels
        bottom_labels = first.col_labels
    else:
        top_labels = first.col_labels
        bottom_labels = None

    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    tick_idx = [ci for ci, _ in top_labels]
    tick_lbl = [lb for _, lb in top_labels]
    ax.set_xticks(tick_idx)
    ax.set_xticklabels(tick_lbl, rotation=0, ha="center", fontproperties=mono_sm)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Right-margin annotation naming the top axis (e.g., "(HxB2)").
    if has_numbering and first.numbering_ref_label:
        ax.text(
            max_cols - 0.2,
            -1.2,
            f"({first.numbering_ref_label})",
            ha="left",
            va="center",
            fontproperties=mono_sm,
            color=PANEL_LABEL_COLOR,
            fontstyle="italic",
        )

    # Bottom secondary x-axis for the variant-call ref.
    if bottom_labels is not None:
        sec = ax.secondary_xaxis("bottom")
        sec_idx = [ci for ci, _ in bottom_labels]
        sec_lbl = [lb for _, lb in bottom_labels]
        sec.set_xticks(sec_idx)
        sec.set_xticklabels(sec_lbl, rotation=0, ha="center", fontproperties=mono_sm)
        sec.tick_params(axis="x", length=0, pad=2)
        for spine in sec.spines.values():
            spine.set_visible(False)
        if first.variant_ref_label:
            ax.text(
                max_cols - 0.2,
                total_rows - 0.5 + 1.2,
                f"({first.variant_ref_label})",
                ha="left",
                va="center",
                fontproperties=mono_sm,
                color=PANEL_LABEL_COLOR,
                fontstyle="italic",
            )

    return ax


def render_panels(
    panels: list[Panel],
    out_path: str | Path = "alignment.png",
    fontsize: float = 12,
    dpi: int = 600,
    palette: str = "nt",
    cell: float | None = None,
    classic: bool = False,
    show_row_labels: bool = False,
    row_label_width: float = 8.0,
    per_panel_axes: bool = False,
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
        per_panel_axes: When True, draw inline tick-label rows per panel
            instead of placing tick labels only at the figure edges.
    """
    fig_w, fig_h = panel_figsize(panels, fontsize, cell, per_panel_axes=per_panel_axes)
    if show_row_labels:
        if cell is None:
            cell = fontsize / 72
        fig_w += row_label_width * cell
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    draw_panels(
        panels,
        ax,
        fontsize=fontsize,
        palette=palette,
        cell=cell,
        classic=classic,
        show_row_labels=show_row_labels,
        row_label_width=row_label_width,
        per_panel_axes=per_panel_axes,
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
    log = logging.getLogger(__name__)
    log.info(
        "Saved: %s (%d dpi, %d panel(s), %d cols)", out_path, dpi, len(panels), max_cols
    )
