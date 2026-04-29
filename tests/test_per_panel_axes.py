"""Tests for the --per-panel-axes flag (inline per-panel x-axis tick rows)."""

from __future__ import annotations

from click.testing import CliRunner

from tview.cli import main
from tview.fasta import fasta_panel
from tview.renderer import panel_figsize, render_panels

from .conftest import OUTPUT_DIR


def _dual_ref_panels(write_fasta, n: int) -> list:
    """Build *n* dual-ref panels (HxB2 numbering + SF162 variant)."""
    panels = []
    for i in range(n):
        path = write_fasta(
            [
                ("HxB2", "ACGT--ACGTACG"),
                ("SF162", "ACGTACG-ACGTA"),
                (f"sample_{i}_a", "ACGTACGTACGTA"),
                (f"sample_{i}_b", "ACGAACGTACGTA"),
            ],
            filename=f"dual_{i}.fasta",
        )
        panels.append(fasta_panel(str(path), variant_ref="SF162", numbering_ref="HxB2"))
    return panels


def _single_ref_panels(write_fasta, n: int) -> list:
    """Build *n* single-ref panels (no numbering layering)."""
    panels = []
    for i in range(n):
        path = write_fasta(
            [
                ("ref", "ACGTACGTAC"),
                (f"s{i}_a", "ACGTACGTAC"),
                (f"s{i}_b", "ACGAACGTAC"),
            ],
            filename=f"single_{i}.fasta",
        )
        panels.append(fasta_panel(str(path)))
    return panels


class TestPerPanelAxesRender:
    def test_single_panel_flag_no_op(self, write_fasta):
        """Single panel + flag on still renders a non-empty file."""
        panel = _dual_ref_panels(write_fasta, 1)[0]
        out = OUTPUT_DIR / "per_panel_axes_single.png"
        render_panels(
            [panel],
            str(out),
            palette="aa",
            dpi=120,
            per_panel_axes=True,
        )
        assert out.exists()
        assert out.stat().st_size > 0

    def test_multi_panel_dual_ref_extra_rows(self, write_fasta):
        """Dual-ref multi-panel: figure height grows when flag is on."""
        panels = _dual_ref_panels(write_fasta, 3)
        w_off, h_off = panel_figsize(panels, fontsize=12)
        w_on, h_on = panel_figsize(panels, fontsize=12, per_panel_axes=True)
        # Each dual-ref panel adds 2 rows (top + bottom tick) when flag on.
        assert w_on == w_off
        assert h_on > h_off

        out = OUTPUT_DIR / "per_panel_axes_multi_dual.png"
        render_panels(
            panels,
            str(out),
            palette="aa",
            dpi=120,
            per_panel_axes=True,
        )
        assert out.exists()
        assert out.stat().st_size > 0

    def test_multi_panel_single_ref_extra_rows(self, write_fasta):
        """Single-ref multi-panel: figure height grows when flag is on (top tick row only)."""
        panels = _single_ref_panels(write_fasta, 3)
        _, h_off = panel_figsize(panels, fontsize=12)
        _, h_on = panel_figsize(panels, fontsize=12, per_panel_axes=True)
        # Each single-ref panel adds 1 row (top tick) when flag on.
        assert h_on > h_off

        out = OUTPUT_DIR / "per_panel_axes_multi_single.png"
        render_panels(
            panels,
            str(out),
            palette="nt",
            dpi=120,
            per_panel_axes=True,
        )
        assert out.exists()
        assert out.stat().st_size > 0


class TestPerPanelAxesVisual:
    def test_multi_panel_dual_ref_visual_fixture(self, write_fasta, output_dir):
        """3-panel dual-ref figure with redundant per-panel x-axes.

        Produces a visual fixture in tests/output/ for manual inspection. Each
        of the three animal panels should show its own HxB2 tick row above and
        SF162p3_ref tick row below, with right-margin axis labels.
        """
        animals = [
            ("RM_A", ["AAGTAGGAGTACAACA", "AAGAACGAGTACAGCA", "AAGTACGAGTACAACA"]),
            ("RM_B", ["AAGTAAGAGTACAACA", "AAGTAGGAGTACATCA", "AAGTAGGTGTACAACA"]),
            ("RM_C", ["AAGCAGGAGTACAACA", "AAGTAGGAGCACAACA", "AAGTGGGAGTACAACA"]),
        ]
        panels = []
        for animal_id, samples in animals:
            seqs = [
                ("HxB2", "AAGTAGGAGTAC--CA"),
                ("SF162p3_ref", "AAGTAGGAGTACA-CA"),
            ]
            for j, s in enumerate(samples):
                seqs.append((f"{animal_id}_clone_{j+1}", s))
            path = write_fasta(seqs, filename=f"{animal_id}.fasta")
            panels.append(
                fasta_panel(
                    str(path),
                    variant_ref="SF162p3_ref",
                    numbering_ref="HxB2",
                )
            )

        out = output_dir / "per_panel_axes_dual_ref_multi.png"
        render_panels(
            panels,
            str(out),
            palette="aa",
            dpi=150,
            show_row_labels=True,
            per_panel_axes=True,
            fontsize=10,
        )
        assert out.exists()
        assert out.stat().st_size > 0

        # Sanity: figure with the flag is taller than without.
        _, h_off = panel_figsize(panels, fontsize=10)
        _, h_on = panel_figsize(panels, fontsize=10, per_panel_axes=True)
        assert h_on > h_off


class TestPerPanelAxesCLI:
    def test_cli_per_panel_axes(self, write_fasta):
        """CLI invocation with --per-panel-axes writes a non-empty output."""
        path = write_fasta(
            [
                ("HxB2", "ACGT--ACGTACG"),
                ("SF162", "ACGTACG-ACGTA"),
                ("CH505_a", "ACGTACGTACGTA"),
            ]
        )
        out = OUTPUT_DIR / "cli_per_panel_axes.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(path),
                "--fasta",
                str(path),
                "--variant-ref",
                "SF162",
                "--numbering-ref",
                "HxB2",
                "--per-panel-axes",
                "-o",
                str(out),
                "--palette",
                "aa",
                "--dpi",
                "120",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert out.stat().st_size > 0
