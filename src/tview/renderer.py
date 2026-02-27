"""Rendering engine for alignment panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

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
        import glob as _gl

        candidates = _gl.glob(
            str(Path(fm.findfont("monospace")).parent / FONT_FALLBACK_FILENAME)
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
    """Render alignment panels to a publication-quality image file.

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
    """
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
