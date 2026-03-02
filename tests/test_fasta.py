"""Tests for FASTA parsing and panel construction using Hypothesis."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from tview.fasta import fasta_panel, read_fasta
from tview.models import Panel
from tview.renderer import render_panels

from .conftest import OUTPUT_DIR

# ── Helpers ───────────────────────────────────────────────────────


def _write_fasta(seqs: list[tuple[str, str]]) -> str:
    """Write sequences to a temp FASTA file, return path.

    Uses NamedTemporaryFile for Hypothesis compatibility (no tmp_path fixture).
    """
    lines = "".join(f">{name}\n{seq}\n" for name, seq in seqs)
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    f.write(lines)
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
        """Column slicing with contiguous range produces correct subset width."""
        path = _write_fasta(seqs)
        try:
            aln_len = len(seqs[0][1])
            col_end = min(10, aln_len)
            cols = list(range(1, col_end + 1))
            panel = fasta_panel(path, columns=cols)
            assert panel.total_cols == len(cols)
            assert len(panel.ref_row) == len(cols)
        finally:
            Path(path).unlink(missing_ok=True)

    def test_discrete_column_selection(self, tmp_path):
        """Discrete column positions select only those columns."""
        path = tmp_path / "disc.fasta"
        path.write_text(">ref\nABCDEFGHIJ\n>s1\nabcdefghij\n")
        panel = fasta_panel(str(path), columns=[1, 5, 10])
        assert panel.total_cols == 3
        assert panel.ref_row == ["A", "E", "J"]
        _, row, _ = panel.seq_rows[0]
        assert row == ["A", "E", "J"]

    def test_mixed_range_and_discrete(self, tmp_path):
        """Mixed ranges and discrete positions produce correct width."""
        path = tmp_path / "mixed.fasta"
        path.write_text(">ref\nABCDEFGHIJKL\n>s1\nabcdefghijkl\n")
        # columns 1-4, 8, 12 → 6 columns total
        cols = [1, 2, 3, 4, 8, 12]
        panel = fasta_panel(str(path), columns=cols)
        assert panel.total_cols == 6
        assert panel.ref_row == ["A", "B", "C", "D", "H", "L"]

    def test_out_of_range_columns_ignored(self, tmp_path):
        """Column positions beyond sequence length are silently ignored."""
        path = tmp_path / "short.fasta"
        path.write_text(">ref\nACGT\n>s1\nACGA\n")
        panel = fasta_panel(str(path), columns=[1, 2, 100])
        assert panel.total_cols == 2
        assert panel.ref_row == ["A", "C"]

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

    def test_classic_mode_nt(self, tmp_path, output_dir):
        """Classic mode renders black-and-white NT alignment without crashing."""
        f = tmp_path / "classic_nt.fasta"
        f.write_text(">ref\nACGTACGT\n>s1\nACCTACGA\n>s2\nACGTACGT\n")
        panel = fasta_panel(str(f))
        out = output_dir / "classic_mode_nt.png"
        render_panels([panel], str(out), palette="nt", dpi=150, classic=True)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_classic_mode_aa(self, tmp_path, output_dir):
        """Classic mode renders black-and-white AA alignment without crashing."""
        f = tmp_path / "classic_aa.fasta"
        f.write_text(">ref\nMRVKEKYQ\n>s1\nMRVGEKYQ\n>s2\nMRVKDKYQ\n")
        panel = fasta_panel(str(f))
        out = output_dir / "classic_mode_aa.png"
        render_panels([panel], str(out), palette="aa", dpi=150, classic=True)
        assert out.exists()
        assert out.stat().st_size > 0


# ── Panel dataclass sanity ────────────────────────────────────────


class TestDrawPanels:
    """Tests for draw_panels() function."""

    def test_draws_on_provided_axes(self, tmp_path):
        """draw_panels returns the same axes object it was given."""
        from tview.renderer import draw_panels

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGT\n>s1\nACGA\n")
        panel = fasta_panel(str(f))
        fig, ax = plt.subplots(figsize=(4, 1))
        result = draw_panels([panel], ax)
        assert result is ax
        plt.close(fig)

    def test_axes_limits_set_correctly(self, tmp_path):
        """draw_panels sets xlim and ylim based on panel dimensions."""
        from tview.renderer import draw_panels

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n>s2\nACGTACGA\n")
        panel = fasta_panel(str(f))
        fig, ax = plt.subplots(figsize=(4, 2))
        draw_panels([panel], ax)
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        assert xlim == (-0.5, 7.5)
        assert ylim[0] > ylim[1]  # y is inverted
        plt.close(fig)

    def test_no_file_created(self, tmp_path):
        """draw_panels does not create any files."""
        from tview.renderer import draw_panels

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGT\n>s1\nACGA\n")
        panel = fasta_panel(str(f))
        fig, ax = plt.subplots(figsize=(4, 1))
        files_before = set(tmp_path.iterdir())
        draw_panels([panel], ax)
        files_after = set(tmp_path.iterdir())
        assert files_before == files_after
        plt.close(fig)

    def test_classic_mode(self, tmp_path):
        """draw_panels with classic=True does not crash."""
        from tview.renderer import draw_panels

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGT\n>s1\nACGA\n")
        panel = fasta_panel(str(f))
        fig, ax = plt.subplots(figsize=(4, 1))
        result = draw_panels([panel], ax, classic=True)
        assert result is ax
        plt.close(fig)

    def test_aa_palette(self, tmp_path):
        """draw_panels with palette='aa' does not crash."""
        from tview.renderer import draw_panels

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nMRVK\n>s1\nMRGK\n")
        panel = fasta_panel(str(f))
        fig, ax = plt.subplots(figsize=(4, 1))
        result = draw_panels([panel], ax, palette="aa")
        assert result is ax
        plt.close(fig)

    def test_stacked_panels(self, tmp_path):
        """draw_panels with multiple panels on a single external axes."""
        from tview.renderer import draw_panels, panel_figsize

        f1 = tmp_path / "g1.fasta"
        f2 = tmp_path / "g2.fasta"
        f1.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n")
        f2.write_text(">ref\nACGTACGT\n>s2\nACGTACGA\n")
        p1 = fasta_panel(str(f1))
        p2 = fasta_panel(str(f2))
        w, h = panel_figsize([p1, p2])
        fig, ax = plt.subplots(figsize=(w, h))
        result = draw_panels([p1, p2], ax)
        assert result is ax
        plt.close(fig)

    @given(aligned_nt_seqs(min_seqs=2, max_seqs=6, min_len=10, max_len=50))
    @settings(max_examples=10, deadline=30000)
    def test_no_crash_hypothesis(self, seqs):
        """Rendering random NT alignments onto external axes does not crash."""
        from tview.renderer import draw_panels, panel_figsize

        fasta_path = _write_fasta(seqs)
        try:
            panel = fasta_panel(fasta_path)
            w, h = panel_figsize([panel])
            fig, ax = plt.subplots(figsize=(w, h))
            result = draw_panels([panel], ax)
            assert result is ax
            plt.close(fig)
        finally:
            Path(fasta_path).unlink(missing_ok=True)

    def test_savefig_after_draw(self, tmp_path, output_dir):
        """draw_panels output can be saved manually by the caller."""
        from tview.renderer import draw_panels, panel_figsize

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n")
        panel = fasta_panel(str(f))
        w, h = panel_figsize([panel])
        fig, ax = plt.subplots(figsize=(w, h))
        draw_panels([panel], ax)
        out = output_dir / "draw_panels_manual_save.png"
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        assert out.exists()
        assert out.stat().st_size > 0


class TestPanelFigsize:
    """Tests for panel_figsize() function."""

    def test_returns_tuple(self, tmp_path):
        """panel_figsize returns a (width, height) tuple."""
        from tview.renderer import panel_figsize

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGT\n>s1\nACGA\n")
        panel = fasta_panel(str(f))
        result = panel_figsize([panel])
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert all(isinstance(v, (int, float)) for v in result)

    def test_minimum_width(self, tmp_path):
        """Figure width is at least 4 inches."""
        from tview.renderer import panel_figsize

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nAC\n>s1\nAC\n")
        panel = fasta_panel(str(f))
        w, _h = panel_figsize([panel])
        assert w >= 4

    def test_minimum_height(self, tmp_path):
        """Figure height is at least 1.0 inches."""
        from tview.renderer import panel_figsize

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nAC\n>s1\nAC\n")
        panel = fasta_panel(str(f))
        _w, h = panel_figsize([panel])
        assert h >= 1.0

    def test_wider_alignment_larger_width(self, tmp_path):
        """Wider alignments produce wider figure sizes."""
        from tview.renderer import panel_figsize

        short = tmp_path / "short.fasta"
        short.write_text(">ref\nACGT\n>s1\nACGA\n")
        long_ = tmp_path / "long.fasta"
        long_.write_text(
            ">ref\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            ">s1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA\n"
        )
        p_short = fasta_panel(str(short))
        p_long = fasta_panel(str(long_))
        w_short, _ = panel_figsize([p_short], cell=0.2)
        w_long, _ = panel_figsize([p_long], cell=0.2)
        assert w_long > w_short

    def test_consistent_with_render_panels(self, tmp_path):
        """panel_figsize produces the same dimensions render_panels would compute."""
        from tview.renderer import panel_figsize

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n")
        panel = fasta_panel(str(f))
        w, h = panel_figsize([panel], fontsize=12, cell=None)
        cell = 12 / 72
        max_cols = panel.total_cols
        total_rows = 1 + len(panel.seq_rows)
        expected_w = max(4, max_cols * cell + 0.5)
        expected_h = max(1.0, total_rows * cell + 0.6)
        assert w == pytest.approx(expected_w)
        assert h == pytest.approx(expected_h)


try:
    import patchworklib  # noqa: F401

    _has_patchworklib = True
except ImportError:
    _has_patchworklib = False


@pytest.mark.skipif(not _has_patchworklib, reason="patchworklib not installed")
class TestPatchworklib:
    """Integration tests using patchworklib Bricks."""

    def test_draw_on_brick(self, tmp_path, output_dir):
        """draw_panels renders onto a patchworklib Brick."""
        import patchworklib as pw

        from tview.renderer import draw_panels, panel_figsize

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n>s2\nACGTACGA\n")
        panel = fasta_panel(str(f))
        w, h = panel_figsize([panel], fontsize=7, cell=0.14)
        brick = pw.Brick(label="alignment", figsize=(w, h))
        result = draw_panels([panel], ax=brick, fontsize=7, cell=0.14)
        assert result is brick
        out = output_dir / "patchworklib_single_brick.png"
        brick.savefig(out, dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_compose_vertical(self, tmp_path, output_dir):
        """Two tview Bricks composed vertically with / operator."""
        import patchworklib as pw

        from tview.renderer import draw_panels, panel_figsize

        f1 = tmp_path / "g1.fasta"
        f2 = tmp_path / "g2.fasta"
        f1.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n")
        f2.write_text(">ref\nACGTACGT\n>s2\nACGTACGA\n")
        p1 = fasta_panel(str(f1))
        p2 = fasta_panel(str(f2))

        w1, h1 = panel_figsize([p1], fontsize=7, cell=0.14)
        brick1 = pw.Brick(label="group1", figsize=(w1, h1))
        draw_panels([p1], ax=brick1, fontsize=7, cell=0.14)

        w2, h2 = panel_figsize([p2], fontsize=7, cell=0.14)
        brick2 = pw.Brick(label="group2", figsize=(w2, h2))
        draw_panels([p2], ax=brick2, fontsize=7, cell=0.14)

        layout = brick1 / brick2
        out = output_dir / "patchworklib_vertical.png"
        layout.savefig(out, dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_compose_with_matplotlib(self, tmp_path, output_dir):
        """tview Brick composed horizontally with a plain matplotlib Brick."""
        import patchworklib as pw

        from tview.renderer import draw_panels, panel_figsize

        f = tmp_path / "test.fasta"
        f.write_text(">ref\nACGTACGT\n>s1\nACCTACGT\n")
        panel = fasta_panel(str(f))
        w, h = panel_figsize([panel], fontsize=7, cell=0.14)

        alignment = pw.Brick(label="aln", figsize=(w, h))
        draw_panels([panel], ax=alignment, fontsize=7, cell=0.14)

        scatter = pw.Brick(label="scatter", figsize=(3, 3))
        scatter.scatter([1, 2, 3], [4, 5, 6])

        layout = alignment | scatter
        out = output_dir / "patchworklib_mixed.png"
        layout.savefig(out, dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0


class TestPanel:
    def test_defaults(self):
        p = Panel("test", ["A", "C"], [], 2, [])
        assert p.ins_columns == set()
        assert p.label == "test"
