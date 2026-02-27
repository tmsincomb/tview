"""Configuration loader for tview color palettes and style settings.

Reads palette.yaml and style.yaml from the package directory and exposes
backward-compatible module-level constants.

Examples:
    >>> palette = load_palette()
    >>> sorted(palette.keys()) == ['amino_acid', 'nucleotide', 'rendering']
    True
    >>> style = load_style()
    >>> sorted(style.keys()) == ['alpha', 'font']
    True
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

_PKG_DIR = Path(__file__).parent


@lru_cache(maxsize=1)
def load_palette() -> dict[str, Any]:
    """Load color palette definitions from palette.yaml.

    Returns:
        Parsed YAML dict with ``nucleotide``, ``amino_acid``, and ``rendering`` keys.

    Examples:
        >>> palette = load_palette()
        >>> 'A' in palette['nucleotide']
        True
        >>> 'mismatch_bg' in palette['rendering']
        True
    """
    return yaml.safe_load((_PKG_DIR / "palette.yaml").read_text())


@lru_cache(maxsize=1)
def load_style() -> dict[str, Any]:
    """Load style settings from style.yaml.

    Returns:
        Parsed YAML dict with ``alpha`` and ``font`` keys.

    Examples:
        >>> style = load_style()
        >>> 0.0 <= style['alpha']['forward'] <= 1.0
        True
        >>> isinstance(style['font']['preferences'], list)
        True
    """
    return yaml.safe_load((_PKG_DIR / "style.yaml").read_text())


# -- Backward-compatible constants -----------------------------------------
_palette = load_palette()
_style = load_style()

NT_COLORS: dict[str, str] = _palette["nucleotide"]
AA_COLORS: dict[str, str] = _palette["amino_acid"]
MISMATCH_BG: str = _palette["rendering"]["mismatch_bg"]
INS_BG: str = _palette["rendering"]["insertion_bg"]
GAP_COLOR: str = _palette["rendering"]["gap_color"]
SEPARATOR_COLOR: str = _palette["rendering"]["separator_color"]
TEXT_COLOR: str = _palette["rendering"]["text_color"]
FALLBACK_BASE_COLOR: str = _palette["rendering"]["fallback_base_color"]
PANEL_LABEL_COLOR: str = _palette["rendering"]["panel_label_color"]

FWD_ALPHA: float = _style["alpha"]["forward"]
REV_ALPHA: float = _style["alpha"]["reverse"]
FONT_PREFERENCES: list[dict[str, str]] = _style["font"]["preferences"]
FONT_FALLBACK_FILENAME: str = _style["font"]["fallback_filename"]
