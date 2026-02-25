"""tview — Publication-quality alignment viewer."""

from importlib.metadata import PackageNotFoundError, version

from tview.tview import (
    AA_COLORS,
    NT_COLORS,
    Panel,
    bam_panel,
    fasta_panel,
    read_fasta,
    render_panels,
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
