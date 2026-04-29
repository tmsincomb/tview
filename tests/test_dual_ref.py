"""Tests for the dual-reference (variant + numbering) feature."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from tview.cli import main
from tview.fasta import fasta_panel
from tview.renderer import render_panels

from .conftest import OUTPUT_DIR

ENV_FASTA = Path(__file__).parent / "data" / "env_protein_aligned.fasta"


# ── fasta_panel: parsing ──────────────────────────────────────────


class TestFastaPanelDualRef:
    def test_variant_ref_default(self, write_fasta):
        """No flags: first seq becomes ref_row, no numbering layering."""
        path = write_fasta([("ref", "ACGT"), ("s1", "ACGA"), ("s2", "ATGT")])
        panel = fasta_panel(str(path))
        assert panel.ref_row == ["A", "C", "G", "T"]
        assert panel.numbering_ref_row is None
        assert panel.numbering_col_labels is None
        assert panel.numbering_ref_label is None
        assert panel.variant_ref_label is None
        assert len(panel.seq_rows) == 2

    def test_variant_ref_by_name(self, write_fasta):
        """variant_ref selects a different sequence as ref_row."""
        path = write_fasta([("ref0", "AAAA"), ("target", "TTTT"), ("sample", "TAAA")])
        panel = fasta_panel(str(path), variant_ref="target")
        assert panel.ref_row == ["T", "T", "T", "T"]
        # Both ref0 and sample remain in seq_rows; target is excluded.
        names = [n for n, _r, _rev in panel.seq_rows]
        assert "target" not in names
        assert "ref0" in names
        assert "sample" in names
        # No numbering layering when only variant_ref is given.
        assert panel.numbering_ref_row is None

    def test_variant_ref_not_found_raises(self, write_fasta):
        """Missing variant_ref name raises ValueError with available headers."""
        path = write_fasta([("ref", "ACGT"), ("s1", "ACGA")])
        with pytest.raises(ValueError, match="variant-ref 'MISSING' not found"):
            fasta_panel(str(path), variant_ref="MISSING")

    def test_numbering_ref_only(self, write_fasta):
        """Only numbering_ref given: first seq stays as variant ref."""
        path = write_fasta([("first", "ACGT"), ("HxB2", "ACGT"), ("sample", "ACGA")])
        panel = fasta_panel(str(path), numbering_ref="HxB2")
        assert panel.ref_row == ["A", "C", "G", "T"]
        assert panel.numbering_ref_row == ["A", "C", "G", "T"]
        assert panel.numbering_ref_label == "HxB2"
        assert panel.variant_ref_label == "first"
        names = [n for n, _r, _rev in panel.seq_rows]
        assert "HxB2" not in names
        assert "first" not in names
        assert "sample" in names

    def test_variant_and_numbering_same_name(self, write_fasta):
        """Both flags equal → collapses to single-ref mode."""
        path = write_fasta([("HxB2", "ACGT"), ("s1", "ACGA")])
        panel = fasta_panel(str(path), variant_ref="HxB2", numbering_ref="HxB2")
        assert panel.ref_row == ["A", "C", "G", "T"]
        assert panel.numbering_ref_row is None
        assert panel.numbering_col_labels is None

    def test_variant_and_numbering_different(self, write_fasta):
        """Both set, different names → both fields populated, neither in seq_rows."""
        # HxB2 has gaps at columns 4-5 (after ACGT); SF162 has gap at 7.
        # → numbering by HxB2 non-gap differs from numbering by SF162 non-gap.
        path = write_fasta(
            [
                ("HxB2", "ACGT--ACGTACG"),
                ("SF162", "ACGTACG-ACGTA"),
                ("CH505", "ACGTACGTACGTA"),
            ]
        )
        panel = fasta_panel(str(path), variant_ref="SF162", numbering_ref="HxB2")
        assert panel.ref_row == list("ACGTACG-ACGTA")
        assert panel.numbering_ref_row == list("ACGT--ACGTACG")
        assert panel.numbering_ref_label == "HxB2"
        assert panel.variant_ref_label == "SF162"
        names = [n for n, _r, _rev in panel.seq_rows]
        assert names == ["CH505"]
        # Different gap patterns → different non-gap position counts at every
        # column. The two label sets should not be identical.
        assert panel.col_labels != panel.numbering_col_labels

    def test_max_rows_with_two_refs(self, write_fasta):
        """max_rows counts samples only when refs are explicitly chosen."""
        seqs = [("HxB2", "ACGT"), ("SF162", "ACGT")] + [
            (f"s{i}", "ACGT") for i in range(10)
        ]
        path = write_fasta(seqs)
        panel = fasta_panel(
            str(path),
            variant_ref="SF162",
            numbering_ref="HxB2",
            max_rows=3,
        )
        assert len(panel.seq_rows) == 3
        names = [n for n, _r, _rev in panel.seq_rows]
        assert names == ["s0", "s1", "s2"]

    def test_max_rows_single_ref_unchanged(self, write_fasta):
        """No new flags + max_rows preserves today's behavior (refs counted)."""
        seqs = [("ref", "ACGT")] + [(f"s{i}", "ACGT") for i in range(10)]
        path = write_fasta(seqs)
        panel = fasta_panel(str(path), max_rows=3)
        assert len(panel.seq_rows) == 3
        assert [n for n, _r, _rev in panel.seq_rows] == ["s0", "s1", "s2"]

    def test_columns_with_numbering_ref(self, write_fasta):
        """--columns subsetting still works in dual-ref mode."""
        path = write_fasta(
            [
                ("HxB2", "ACGT--ACGT"),
                ("SF162", "ACGTACGTAC"),
                ("CH505", "ACGTACGTAC"),
            ]
        )
        panel = fasta_panel(
            str(path),
            columns=[1, 2, 3, 4, 5],
            variant_ref="SF162",
            numbering_ref="HxB2",
        )
        assert panel.total_cols == 5
        assert panel.ref_row == ["A", "C", "G", "T", "A"]
        assert panel.numbering_ref_row == ["A", "C", "G", "T", "-"]
        # numbering_col_labels skip the gap column.
        ncl = panel.numbering_col_labels
        assert ncl is not None
        # First non-gap → label "1" at column 0; gap at column 4 → no label there.
        idx_with_labels = {ci for ci, _ in ncl}
        assert 4 not in idx_with_labels


# ── Renderer smoke ────────────────────────────────────────────────


class TestDualRefRender:
    def test_dual_ref_render_smoke(self, write_fasta, output_dir):
        """End-to-end: build panel with dual refs and render to PNG."""
        path = write_fasta(
            [
                ("HxB2", "ACGT--ACGTACG"),
                ("SF162", "ACGTACG-ACGTA"),
                ("CH505_a", "ACGTACGTACGTA"),
                ("CH505_b", "ACGAACGTACGTA"),
            ]
        )
        panel = fasta_panel(str(path), variant_ref="SF162", numbering_ref="HxB2")
        out = output_dir / "dual_ref_synthetic.png"
        render_panels([panel], str(out), palette="aa", dpi=120, show_row_labels=True)
        assert out.exists()
        assert out.stat().st_size > 0

    @pytest.mark.skipif(
        not ENV_FASTA.exists(), reason="env_protein_aligned.fasta not present"
    )
    def test_dual_ref_render_env_protein(self, output_dir):
        """Realistic HIV Env: HxB2 numbering + SF162p3_ref variant calls."""
        panel = fasta_panel(
            str(ENV_FASTA),
            columns=list(range(1, 61)),
            variant_ref="SF162p3_ref",
            numbering_ref="HxB2",
            max_rows=4,
        )
        assert panel.numbering_ref_row is not None
        assert panel.variant_ref_label == "SF162p3_ref"
        assert panel.numbering_ref_label == "HxB2"
        out = output_dir / "dual_ref_HxB2_SF162p3.png"
        render_panels(
            [panel],
            str(out),
            palette="aa",
            dpi=150,
            show_row_labels=True,
            fontsize=6,
            cell=0.12,
        )
        assert out.exists()
        assert out.stat().st_size > 0


# ── CLI ───────────────────────────────────────────────────────────


class TestDualRefCLI:
    def test_cli_variant_ref_missing_errors(self, write_fasta):
        path = write_fasta([("ref", "ACGT"), ("s1", "ACGA")])
        out = OUTPUT_DIR / "cli_dual_missing.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(path),
                "--variant-ref",
                "MISSING",
                "-o",
                str(out),
                "--dpi",
                "100",
            ],
        )
        assert result.exit_code != 0
        assert "MISSING" in result.output

    def test_cli_dual_ref_render(self, write_fasta):
        path = write_fasta(
            [
                ("HxB2", "ACGT--ACGTACG"),
                ("SF162", "ACGTACG-ACGTA"),
                ("CH505_a", "ACGTACGTACGTA"),
            ]
        )
        out = OUTPUT_DIR / "cli_dual_ref.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(path),
                "--variant-ref",
                "SF162",
                "--numbering-ref",
                "HxB2",
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

    def test_cli_show_row_labels_flag(self, write_fasta):
        path = write_fasta([("ref", "ACGT"), ("s1", "ACGA")])
        out = OUTPUT_DIR / "cli_show_row_labels.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(path),
                "--show-row-labels",
                "-o",
                str(out),
                "--dpi",
                "100",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
