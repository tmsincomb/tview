"""Click CLI for tview."""

from __future__ import annotations

import sys

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
def main(
    bam, ref, region, fasta, columns, output, palette, dpi, fontsize, cell, classic_mode
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
        for bam_path in bam_paths:
            panels.append(bam_panel(bam_path, ref, region))

    if fasta_paths:
        cols = parse_columns(columns) if columns else None
        for fasta_path in fasta_paths:
            panels.append(fasta_panel(fasta_path, columns=cols))

    render_panels(
        panels,
        output,
        fontsize=fontsize,
        dpi=dpi,
        palette=palette,
        cell=cell,
        classic=classic_mode,
    )
