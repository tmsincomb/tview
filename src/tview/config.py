"""Configuration loader for tview color palettes and style settings.

Reads palette.yaml and style.yaml from the package directory and exposes
backward-compatible module-level constants.
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
    """
    with open(_PKG_DIR / "palette.yaml") as fh:
        return yaml.safe_load(fh)


@lru_cache(maxsize=1)
def load_style() -> dict[str, Any]:
    """Load style settings from style.yaml.

    Returns:
        Parsed YAML dict with ``alpha`` and ``font`` keys.
    """
    with open(_PKG_DIR / "style.yaml") as fh:
        return yaml.safe_load(fh)


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
