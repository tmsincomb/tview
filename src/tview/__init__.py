"""tview-image — Publication-quality alignment viewer."""

from tview.tview import (
    Panel,
    fasta_panel,
    bam_panel,
    read_fasta,
    render_panels,
    NT_COLORS,
    AA_COLORS,
)

__version__ = "0.1.0"
__all__ = [
    "Panel",
    "fasta_panel",
    "bam_panel",
    "read_fasta",
    "render_panels",
    "NT_COLORS",
    "AA_COLORS",
]
