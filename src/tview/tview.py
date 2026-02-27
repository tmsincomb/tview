"""Backward-compatible re-export shim.

All logic has moved to focused modules: models, fasta, bam, renderer.
This file preserves ``from tview.tview import ...`` for existing code.

Examples:
    >>> from tview.tview import Panel, read_fasta
"""

from tview.bam import bam_panel, build_read_row  # noqa: F401
from tview.config import AA_COLORS, NT_COLORS  # noqa: F401
from tview.fasta import fasta_panel, read_fasta  # noqa: F401
from tview.models import Panel  # noqa: F401
from tview.renderer import draw_panels, panel_figsize, render_panels  # noqa: F401
