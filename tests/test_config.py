"""Tests for YAML config loading and backward-compatible constants."""

from __future__ import annotations

import re

import pytest

from tview.config import (
    AA_COLORS,
    FALLBACK_BASE_COLOR,
    FONT_FALLBACK_FILENAME,
    FONT_PREFERENCES,
    FWD_ALPHA,
    GAP_COLOR,
    INS_BG,
    MISMATCH_BG,
    NT_COLORS,
    PANEL_LABEL_COLOR,
    REV_ALPHA,
    SEPARATOR_COLOR,
    TEXT_COLOR,
    load_palette,
    load_style,
)

HEX_RE = re.compile(r"^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$")


class TestLoadPalette:
    def test_has_required_keys(self):
        palette = load_palette()
        assert "nucleotide" in palette
        assert "amino_acid" in palette
        assert "rendering" in palette

    def test_nucleotide_keys(self):
        palette = load_palette()
        nt = palette["nucleotide"]
        for base in ("A", "C", "G", "T", "N", "-"):
            assert base in nt, f"Missing nucleotide key: {base}"

    def test_amino_acid_keys(self):
        palette = load_palette()
        aa = palette["amino_acid"]
        expected = set("AVLIMFWPKRHDESTNGCYX*-")
        assert expected.issubset(set(aa.keys()))

    def test_rendering_keys(self):
        palette = load_palette()
        rendering = palette["rendering"]
        for key in (
            "mismatch_bg",
            "insertion_bg",
            "gap_color",
            "separator_color",
            "text_color",
            "fallback_base_color",
            "panel_label_color",
        ):
            assert key in rendering, f"Missing rendering key: {key}"

    def test_all_colors_are_hex(self):
        palette = load_palette()
        for section in ("nucleotide", "amino_acid"):
            for key, val in palette[section].items():
                assert HEX_RE.match(val), f"{section}.{key} not valid hex: {val}"
        for key, val in palette["rendering"].items():
            assert HEX_RE.match(val), f"rendering.{key} not valid hex: {val}"


class TestLoadStyle:
    def test_has_required_keys(self):
        style = load_style()
        assert "alpha" in style
        assert "font" in style

    def test_alpha_values_are_floats(self):
        style = load_style()
        assert isinstance(style["alpha"]["forward"], float)
        assert isinstance(style["alpha"]["reverse"], float)

    def test_alpha_range(self):
        style = load_style()
        for key in ("forward", "reverse"):
            val = style["alpha"][key]
            assert 0.0 <= val <= 1.0, f"alpha.{key} out of range: {val}"

    def test_font_preferences_structure(self):
        style = load_style()
        prefs = style["font"]["preferences"]
        assert isinstance(prefs, list)
        assert len(prefs) >= 1
        for pref in prefs:
            assert "family" in pref

    def test_fallback_filename(self):
        style = load_style()
        assert style["font"]["fallback_filename"].endswith(".ttf")


class TestBackwardCompatibleConstants:
    def test_nt_colors_values(self):
        assert NT_COLORS["A"] == "#4CAF50"
        assert NT_COLORS["C"] == "#2196F3"
        assert NT_COLORS["G"] == "#FF9800"
        assert NT_COLORS["T"] == "#F44336"
        assert NT_COLORS["N"] == "#9E9E9E"
        assert NT_COLORS["-"] == "#9E9E9E"

    def test_aa_colors_values(self):
        assert AA_COLORS["K"] == "#F44336"
        assert AA_COLORS["D"] == "#E040FB"
        assert AA_COLORS["S"] == "#4CAF50"
        assert AA_COLORS["G"] == "#FF9800"
        assert AA_COLORS["*"] == "#9E9E9E"

    def test_rendering_constants(self):
        assert MISMATCH_BG == "#FFEB3B55"
        assert INS_BG == "#CE93D833"
        assert GAP_COLOR == "#9E9E9E"
        assert SEPARATOR_COLOR == "#BDBDBD"
        assert TEXT_COLOR == "#000000"
        assert FALLBACK_BASE_COLOR == "#9E9E9E"
        assert PANEL_LABEL_COLOR == "#616161"

    def test_alpha_constants(self):
        assert FWD_ALPHA == 1.0
        assert REV_ALPHA == 0.85

    def test_font_constants(self):
        assert isinstance(FONT_PREFERENCES, list)
        assert FONT_PREFERENCES[0]["family"] == "Helvetica"
        assert FONT_FALLBACK_FILENAME == "DejaVuSansMono-Bold.ttf"
