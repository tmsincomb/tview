"""Tests for FASTA parsing and panel construction using Hypothesis."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from tview.fasta import fasta_panel, read_fasta
from tview.models import Panel
from tview.renderer import render_panels

from .conftest import OUTPUT_DIR

# ── Helpers ───────────────────────────────────────────────────────


def _write_fasta(seqs: list[tuple[str, str]]) -> str:
    """Write sequences to a temp FASTA file, return path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    for name, seq in seqs:
        f.write(f">{name}\n{seq}\n")
    f.close()
    return f.name


def _keep(src: str, name: str) -> None:
    """Copy rendered image into tests/output/ for visual inspection."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    shutil.copy2(src, OUTPUT_DIR / name)


# ── Strategies ────────────────────────────────────────────────────

seq_name = st.text(
    alphabet="abcdefghijklmnopqrstuvwxyz0123456789_",
    min_size=1,
    max_size=20,
)


@st.composite
def aligned_nt_seqs(draw, min_seqs=2, max_seqs=10, min_len=5, max_len=200):
    """Generate aligned NT sequences (same length, with optional gaps)."""
    length = draw(st.integers(min_value=min_len, max_value=max_len))
    n_seqs = draw(st.integers(min_value=min_seqs, max_value=max_seqs))
    seqs = []
    for i in range(n_seqs):
        name = draw(seq_name)
        seq = draw(st.text(alphabet="ACGT-", min_size=length, max_size=length))
        seqs.append((f"{name}_{i}", seq))
    return seqs


@st.composite
def aligned_aa_seqs(draw, min_seqs=2, max_seqs=10, min_len=5, max_len=200):
    """Generate aligned AA sequences (same length, with optional gaps)."""
    length = draw(st.integers(min_value=min_len, max_value=max_len))
    n_seqs = draw(st.integers(min_value=min_seqs, max_value=max_seqs))
    seqs = []
    for i in range(n_seqs):
        name = draw(seq_name)
        seq = draw(
            st.text(alphabet="ACDEFGHIKLMNPQRSTVWY-", min_size=length, max_size=length)
        )
        seqs.append((f"{name}_{i}", seq))
    return seqs


# ── read_fasta tests ──────────────────────────────────────────────


class TestReadFasta:
    @given(aligned_nt_seqs(min_seqs=1, max_seqs=8, min_len=5, max_len=100))
    @settings(max_examples=20)
    def test_roundtrip(self, seqs):
        """Writing then reading FASTA recovers all sequences."""
        path = _write_fasta(seqs)
        try:
            result = read_fasta(path)
            assert len(result) == len(seqs)
            for (exp_name, exp_seq), (got_name, got_seq) in zip(seqs, result):
                assert got_name == exp_name
                assert got_seq == exp_seq
        finally:
            Path(path).unlink(missing_ok=True)

    def test_empty_file(self, tmp_path):
        path = tmp_path / "empty.fasta"
        path.write_text("")
        assert read_fasta(str(path)) == []

    def test_single_sequence(self, tmp_path):
        path = tmp_path / "single.fasta"
        path.write_text(">ref\nACGTACGT\n")
        result = read_fasta(str(path))
        assert len(result) == 1
        assert result[0] == ("ref", "ACGTACGT")


# ── fasta_panel tests ─────────────────────────────────────────────


class TestFastaPanel:
    @given(aligned_nt_seqs(min_seqs=2, max_seqs=8, min_len=10, max_len=100))
    @settings(max_examples=20)
    def test_panel_shape(self, seqs):
        """Panel has correct dimensions from random aligned sequences."""
        path = _write_fasta(seqs)
        try:
            panel = fasta_panel(path)
            aln_len = len(seqs[0][1])
            assert panel.total_cols == aln_len
            assert len(panel.ref_row) == aln_len
            assert len(panel.seq_rows) == len(seqs) - 1
        finally:
            Path(path).unlink(missing_ok=True)

    @given(aligned_nt_seqs(min_seqs=2, max_seqs=5, min_len=20, max_len=80))
    @settings(max_examples=15)
    def test_column_slicing(self, seqs):
        """Column slicing produces correct subset width."""
        path = _write_fasta(seqs)
        try:
            aln_len = len(seqs[0][1])
            col_start, col_end = 1, min(10, aln_len)
            panel = fasta_panel(path, col_start=col_start, col_end=col_end)
            expected_len = col_end - col_start + 1
            assert panel.total_cols == expected_len
            assert len(panel.ref_row) == expected_len
        finally:
            Path(path).unlink(missing_ok=True)

    def test_mismatch_detection(self, tmp_path):
        """Mismatches are preserved in seq_rows."""
        path = tmp_path / "mm.fasta"
        path.write_text(">ref\nAAAA\n>s1\nACGA\n")
        panel = fasta_panel(str(path))
        _, row, _ = panel.seq_rows[0]
        assert row[0] == "A"
        assert row[1] == "C"
        assert row[2] == "G"
        assert row[3] == "A"

    def test_gap_handling(self, tmp_path):
        """Gaps in sequences are preserved."""
        path = tmp_path / "gap.fasta"
        path.write_text(">ref\nACGT\n>s1\nA--T\n")
        panel = fasta_panel(str(path))
        _, row, _ = panel.seq_rows[0]
        assert row == ["A", "-", "-", "T"]

    def test_col_labels_skip_gaps(self, tmp_path):
        """Column labels count only non-gap reference positions."""
        path = tmp_path / "gaps_ref.fasta"
        path.write_text(">ref\nA-CGT-ACGT\n>s1\nAACGTAACGT\n")
        panel = fasta_panel(str(path))
        label_positions = {lbl for _, lbl in panel.col_labels}
        assert "1" in label_positions


# ── render_panels tests ───────────────────────────────────────────


class TestRenderPanels:
    _nt_counter = 0
    _aa_counter = 0

    @given(aligned_nt_seqs(min_seqs=2, max_seqs=6, min_len=10, max_len=50))
    @settings(max_examples=10, deadline=30000)
    def test_render_nt_no_crash(self, seqs):
        """Rendering random NT alignments produces a file; kept in tests/output/."""
        fasta_path = _write_fasta(seqs)
        out = fasta_path + ".out.png"
        try:
            panel = fasta_panel(fasta_path)
            render_panels([panel], out, palette="nt", dpi=150)
            assert Path(out).stat().st_size > 0
            TestRenderPanels._nt_counter += 1
            _keep(out, f"hypothesis_nt_{TestRenderPanels._nt_counter:03d}.png")
        finally:
            Path(fasta_path).unlink(missing_ok=True)
            Path(out).unlink(missing_ok=True)

    @given(aligned_aa_seqs(min_seqs=2, max_seqs=6, min_len=10, max_len=50))
    @settings(max_examples=10, deadline=30000)
    def test_render_aa_no_crash(self, seqs):
        """Rendering random AA alignments produces a file; kept in tests/output/."""
        fasta_path = _write_fasta(seqs)
        out = fasta_path + ".out.png"
        try:
            panel = fasta_panel(fasta_path)
            render_panels([panel], out, palette="aa", dpi=150)
            assert Path(out).stat().st_size > 0
            TestRenderPanels._aa_counter += 1
            _keep(out, f"hypothesis_aa_{TestRenderPanels._aa_counter:03d}.png")
        finally:
            Path(fasta_path).unlink(missing_ok=True)
            Path(out).unlink(missing_ok=True)

    def test_stacked_panels(self, tmp_path, output_dir):
        """Stacking multiple panels produces output."""
        f1 = tmp_path / "g1.fasta"
        f2 = tmp_path / "g2.fasta"
        f1.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n")
        f2.write_text(">ref\nACGTACGT\n>s2\nACGTACGA\n")
        p1 = fasta_panel(str(f1))
        p2 = fasta_panel(str(f2))
        out = output_dir / "stacked_panels.png"
        render_panels([p1, p2], str(out), dpi=150)
        assert out.exists()

    def test_output_formats(self, tmp_path, output_dir):
        """SVG and PDF output work."""
        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGT\n>s1\nACGA\n")
        panel = fasta_panel(str(f))
        for ext in ["svg", "pdf"]:
            out = output_dir / f"format_test.{ext}"
            render_panels([panel], str(out), dpi=150)
            assert out.exists()
            assert out.stat().st_size > 0

    def test_all_matches_dot_rendering(self, tmp_path, output_dir):
        """Identical sequences — all dots."""
        f = tmp_path / "match.fasta"
        f.write_text(">ref\nACGTACGTACGT\n>s1\nACGTACGTACGT\n>s2\nACGTACGTACGT\n")
        panel = fasta_panel(str(f))
        out = output_dir / "all_matches.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()

    def test_all_gaps_sequence(self, tmp_path, output_dir):
        """A sequence that is all gaps."""
        f = tmp_path / "allgap.fasta"
        f.write_text(">ref\nACGT\n>s1\n----\n")
        panel = fasta_panel(str(f))
        out = output_dir / "all_gaps.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()

    def test_hiv_env_realistic(self, tmp_path, output_dir):
        """Realistic HIV Env protein alignment snippet."""
        f = tmp_path / "hiv_env.fasta"
        f.write_text(
            ">HxB2\n"
            "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCAS\n"
            ">CAP256_SU\n"
            "MRVKGIQKNWQHLWRWGTLLLGMLMICSATDKLWVTVYYGVPVWKDADTTLFCAS\n"
            ">CH505_TF\n"
            "MRVMGIQRNCQHLWRWGTLILGMLMICSAADKLWVTV-YGVPVWKEAKTTLFCAS\n"
            ">BG505_SOSIP\n"
            "MRVKEKYQHLWRWGWRWGTMLLG--MICSATEKLWVTVY-GVPVWKEATTTLFCAS\n"
        )
        panel = fasta_panel(str(f))
        out = output_dir / "hiv_env_realistic.png"
        render_panels([panel], str(out), palette="aa", dpi=150, fontsize=6, cell=0.12)
        assert out.exists()

    def test_dense_nt_alignment(self, tmp_path, output_dir):
        """Dense NT alignment with many mismatches."""
        import random

        random.seed(42)
        ref = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        lines = [">ref\n" + ref + "\n"]
        for i in range(8):
            mutant = list(ref)
            for j in random.sample(range(len(ref)), k=10):
                mutant[j] = random.choice([b for b in "ACGT" if b != ref[j]])
            lines.append(f">sample_{i:03d}\n{''.join(mutant)}\n")
        f = tmp_path / "dense_nt.fasta"
        f.write_text("".join(lines))
        panel = fasta_panel(str(f))
        out = output_dir / "dense_nt_mismatches.png"
        render_panels([panel], str(out), palette="nt", dpi=150)
        assert out.exists()


# ── Panel dataclass sanity ────────────────────────────────────────


class TestPanel:
    def test_defaults(self):
        p = Panel("test", ["A", "C"], [], 2, [], ins_columns=None)
        assert p.ins_columns == set()
        assert p.label == "test"
