"""tview — Publication-quality alignment viewer.

Examples:
    >>> import tview
    >>> hasattr(tview, '__version__')
    True
"""

from importlib.metadata import PackageNotFoundError, version

from tview.bam import bam_panel
from tview.config import AA_COLORS, NT_COLORS
from tview.fasta import fasta_panel, read_fasta
from tview.models import Panel
from tview.renderer import render_panels

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
