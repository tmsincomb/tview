"""Click CLI for tview."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import click

from tview.bam import bam_panel
from tview.fasta import fasta_panel
from tview.renderer import render_panels


def parse_columns(columns: str) -> list[int]:
    """Parse column spec into sorted 1-based positions.

    Supports individual positions, ranges, and mixed:
      '5,40,690'    → [5, 40, 690]
      '1-120'       → [1, 2, ..., 120]
      '5,10-20,40'  → [5, 10, 11, ..., 20, 40]

    Args:
        columns: Comma-separated column spec string.

    Returns:
        Sorted, deduplicated list of 1-based column positions.

    Examples:
        >>> parse_columns('5,40,690')
        [5, 40, 690]
        >>> parse_columns('1-5')
        [1, 2, 3, 4, 5]
        >>> parse_columns('5,10-12,40')
        [5, 10, 11, 12, 40]
    """
    positions: set[int] = set()
    for part in columns.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-", 1)
            positions.update(range(int(start), int(end) + 1))
        else:
            positions.add(int(part))
    return sorted(positions)


def _show_in_terminal(output_path: str | Path) -> None:
    """Display a rendered image inline via ``kitten icat``.

    Works in terminals that implement the kitty graphics protocol
    (Kitty, Ghostty). Prints a soft warning to stderr if ``kitten`` is
    not on PATH — the file has already been written, so the render
    itself is not considered failed.

    Args:
        output_path: Path to the image file just rendered.
    """
    kitten = shutil.which("kitten")
    if not kitten:
        click.echo(
            f"--show requires 'kitten' on PATH (install kitty). "
            f"Image saved to: {output_path}",
            err=True,
        )
        return
    subprocess.run(
        [kitten, "icat", "--align", "left", str(output_path)],
        check=False,
    )


def _expand_stdin(paths: list[str]) -> list[str]:
    """If paths is ['-'], read file paths from stdin (one per line).

    Args:
        paths: List of file paths or ['-'] to read from stdin.

    Returns:
        Expanded list of file paths.

    Examples:
        >>> _expand_stdin(['file1.fa', 'file2.fa'])
        ['file1.fa', 'file2.fa']
        >>> _expand_stdin([])
        []
    """
    if paths and len(paths) == 1 and paths[0] == "-":
        return [line.strip() for line in sys.stdin if line.strip()]
    return list(paths)


@click.command(
    context_settings={"help_option_names": ["-h", "--help"]},
    epilog="Use '-' to read file paths from stdin, e.g.:\n\n"
    "  find . -name '*.fasta' | tview --fasta - --palette aa -o out.png",
)
@click.option(
    "--bam",
    multiple=True,
    help="BAM file(s) — each becomes a panel. Use '-' for stdin.",
)
@click.option(
    "--ref",
    type=click.Path(exists=True),
    help="Reference FASTA (required for BAM mode).",
)
@click.option("--region", help="Genomic region chr:start-end (required for BAM mode).")
@click.option(
    "--fasta",
    multiple=True,
    help="Aligned FASTA file(s) — each becomes a panel. Use '-' for stdin.",
)
@click.option(
    "--columns",
    help="Column positions for FASTA, 1-based (e.g. 1-120, 5,40,690, or 5,10-20,40).",
)
@click.option(
    "-o",
    "--output",
    default="alignment.png",
    show_default=True,
    help="Output image path.",
)
@click.option(
    "--palette",
    type=click.Choice(["nt", "aa"]),
    default="nt",
    show_default=True,
    help="Color palette.",
)
@click.option(
    "--dpi", type=int, default=300, show_default=True, help="Image resolution."
)
@click.option(
    "--fontsize",
    type=int,
    default=7,
    show_default=True,
    help="Base font size in points.",
)
@click.option(
    "--cell", type=float, default=0.14, show_default=True, help="Cell size in inches."
)
@click.option(
    "--classic-mode",
    is_flag=True,
    default=False,
    help="Black-and-white rendering with no color highlighting.",
)
@click.option(
    "--show",
    is_flag=True,
    default=False,
    help="Display rendered image inline via 'kitten icat' (Kitty/Ghostty).",
)
@click.option(
    "--max-rows",
    type=int,
    default=None,
    help="Cap non-reference rows (FASTA) or reads (BAM) per panel.",
)
@click.option(
    "--variant-ref",
    default=None,
    help="FASTA header name for variant calling (mismatches drawn against this). "
    "Default: first sequence.",
)
@click.option(
    "--numbering-ref",
    default=None,
    help="FASTA header name for x-axis numbering. When different from "
    "--variant-ref, both refs are shown as rows and two x-axes are drawn.",
)
@click.option(
    "--show-row-labels",
    is_flag=True,
    default=False,
    help="Show sequence ID labels on the left of each row.",
)
@click.option(
    "--tick-every",
    type=int,
    default=10,
    show_default=True,
    help="Label every Nth x-axis position. Use 1 to label every column.",
)
@click.option(
    "--per-panel-axes",
    is_flag=True,
    default=False,
    help="In multi-panel mode, draw x-axis tick labels above (and below "
    "for dual-ref panels) each panel instead of only at the figure edges.",
)
def main(
    bam,
    ref,
    region,
    fasta,
    columns,
    output,
    palette,
    dpi,
    fontsize,
    cell,
    classic_mode,
    show,
    max_rows,
    variant_ref,
    numbering_ref,
    show_row_labels,
    tick_every,
    per_panel_axes,
):
    """Publication-quality alignment viewer (BAM or FASTA).

    Supports BAM files (with reference FASTA), pre-aligned FASTA (e.g. MAFFT
    output), and stacking multiple inputs into a single figure.
    """
    bam_paths = _expand_stdin(list(bam))
    fasta_paths = _expand_stdin(list(fasta))

    if not bam_paths and not fasta_paths:
        raise click.UsageError("Provide --bam and/or --fasta")

    panels = []

    if bam_paths:
        if not ref or not region:
            raise click.UsageError("--ref and --region are required for BAM input")
        if variant_ref or numbering_ref:
            click.echo(
                "warning: --variant-ref/--numbering-ref are FASTA-only and were ignored",
                err=True,
            )
        for bam_path in bam_paths:
            panels.append(
                bam_panel(
                    bam_path,
                    ref,
                    region,
                    max_rows=max_rows,
                    tick_every=tick_every,
                )
            )

    if fasta_paths:
        cols = parse_columns(columns) if columns else None
        for fasta_path in fasta_paths:
            try:
                panels.append(
                    fasta_panel(
                        fasta_path,
                        columns=cols,
                        max_rows=max_rows,
                        variant_ref=variant_ref,
                        numbering_ref=numbering_ref,
                        tick_every=tick_every,
                    )
                )
            except ValueError as exc:
                raise click.UsageError(str(exc)) from exc

    # Auto-enable per-row labels when dual-ref mode is active (different refs).
    effective_show_row_labels = show_row_labels or bool(
        variant_ref and numbering_ref and variant_ref != numbering_ref
    )

    render_panels(
        panels,
        output,
        fontsize=fontsize,
        dpi=dpi,
        palette=palette,
        cell=cell,
        classic=classic_mode,
        show_row_labels=effective_show_row_labels,
        per_panel_axes=per_panel_axes,
    )

    if show:
        _show_in_terminal(output)
