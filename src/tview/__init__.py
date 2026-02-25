"""tview — Publication-quality alignment viewer."""

from importlib.metadata import version, PackageNotFoundError

from tview.tview import (
    Panel,
    fasta_panel,
    bam_panel,
    read_fasta,
    render_panels,
    NT_COLORS,
    AA_COLORS,
)

try:
    __version__ = version("tview")
except PackageNotFoundError:
    __version__ = "0.0.0"
__all__ = [
    "Panel",
    "fasta_panel",
    "bam_panel",
    "read_fasta",
    "render_panels",
    "NT_COLORS",
    "AA_COLORS",
]
