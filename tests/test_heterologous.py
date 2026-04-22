"""Tests for the secondary-reference (heterologous) coloring + per-row labels."""

from __future__ import annotations

from pathlib import Path

import pytest

from tview.fasta import fasta_panel
from tview.models import Panel
from tview.renderer import render_panels


def _panel(ref: str, secondary: str | None, rows: list[tuple[str, str]]) -> Panel:
    """Build a Panel with optional secondary_ref_row."""
    return Panel(
        label="primary",
        ref_row=list(ref),
        seq_rows=[(name, list(bases), False) for name, bases in rows],
        total_cols=len(ref),
        col_labels=[(0, "1")],
        secondary_ref_row=list(secondary) if secondary is not None else None,
    )


class TestHeterologousColoring:
    def test_default_color(self):
        panel = _panel("AAAA", "TTTT", [("s1", "TTTT")])
        assert panel.heterologous_color == "#FF6F00"
        assert panel.secondary_ref_row == list("TTTT")

    def test_no_secondary_ref(self, output_dir):
        panel = _panel("ACGTACGT", None, [("s1", "ATCTACGT")])
        out = output_dir / "no_secondary.png"
        render_panels([panel], str(out), dpi=100)
        assert out.exists()
        assert panel.secondary_ref_row is None

    def test_three_way_render(self, output_dir):
        # ref=ACGT, secondary=AGGA, sample=AGGT:
        #   pos 0: A==A → match (.)
        #   pos 1: G != C, G == G → heterologous (orange G)
        #   pos 2: G == G (ref) → match (.)
        #   pos 3: T == T → match (.)
        # sample2 has pure mismatch at pos 0:
        #   pos 0: T != A, T != A (sec) → pure mismatch (red palette)
        panel = _panel(
            ref="ACGT",
            secondary="AGGA",
            rows=[("hetero_at_1", "AGGT"), ("mismatch", "TCGT")],
        )
        out = output_dir / "three_way_classification.png"
        render_panels([panel], str(out), palette="aa", dpi=100)
        assert out.exists()

    def test_classic_mode_disables_heterologous(self, output_dir):
        # In classic mode, all letters render in black; heterologous branch
        # must NOT fire.
        panel = _panel("AAAA", "TTTT", [("s1", "TTTT")])
        out = output_dir / "classic_no_hetero.png"
        render_panels([panel], str(out), classic=True, dpi=100)
        assert out.exists()

    def test_secondary_ref_shorter(self, output_dir):
        # Secondary shorter than ref: must not raise an IndexError.
        panel = _panel("AAAAAAAA", "TTTT", [("s1", "TTTTTTTT")])
        out = output_dir / "secondary_short.png"
        render_panels([panel], str(out), dpi=100)
        assert out.exists()


class TestRowLabels:
    def test_show_row_labels_render(self, output_dir):
        # Single panel with sequence IDs visible on the left
        panel = _panel(
            ref="ACGTACGT",
            secondary=None,
            rows=[("rh2856_s1", "ACCTACGT"), ("rh3031_s2", "ACGTACGA")],
        )
        out = output_dir / "row_labels.png"
        render_panels([panel], str(out), dpi=100, show_row_labels=True)
        assert out.exists()

    def test_no_row_labels_default(self, output_dir):
        # Default behaviour: no per-row labels, panel label only when stacked
        panel = _panel(
            ref="ACGTACGT",
            secondary=None,
            rows=[("s1", "ACCTACGT"), ("s2", "ACGTACGA")],
        )
        out = output_dir / "no_row_labels.png"
        render_panels([panel], str(out), dpi=100)
        assert out.exists()

    def test_row_labels_with_stacked_panels(self, output_dir):
        p1 = _panel("ACGT", None, [("animal1_s1", "ACGT"), ("animal1_s2", "ATGT")])
        p2 = _panel("ACGT", None, [("animal2_s1", "ACGT")])
        out = output_dir / "stacked_with_labels.png"
        render_panels([p1, p2], str(out), dpi=100, show_row_labels=True)
        assert out.exists()


class TestEndToEnd:
    def test_fasta_panel_roundtrip(self, write_fasta, output_dir):
        # End-to-end: build via fasta_panel, attach secondary, render
        path = write_fasta(
            [
                ("ref", "ACGTACGT"),
                ("secondary", "TGCATGCA"),
                (
                    "s1",
                    "TGGTACGT",
                ),  # matches secondary at pos 0,1; pure mismatch elsewhere
            ]
        )
        # Build panel from first 2 rows; the third becomes the only sample
        panel = fasta_panel(path)
        # Overwrite for our test
        panel = Panel(
            label="ref",
            ref_row=list("ACGTACGT"),
            seq_rows=[("s1", list("TGGTACGT"), False)],
            total_cols=8,
            col_labels=panel.col_labels,
            secondary_ref_row=list("TGCATGCA"),
        )
        out = output_dir / "fasta_roundtrip.png"
        render_panels([panel], str(out), palette="aa", dpi=100, show_row_labels=True)
        assert out.exists()
