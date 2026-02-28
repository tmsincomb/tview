"""Tests for BAM parsing and panel construction."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from tview.bam import (
    CIGAR_DEL,
    CIGAR_INS,
    CIGAR_MATCH,
    CIGAR_REF_SKIP,
    CIGAR_SEQ_MATCH,
    CIGAR_SEQ_MISMATCH,
    CIGAR_SOFT_CLIP,
    build_read_row,
)

# ── Helpers ───────────────────────────────────────────────────────


def _mock_read(
    query_name: str,
    reference_start: int,
    query_sequence: str,
    cigartuples: list[tuple[int, int]],
    is_reverse: bool = False,
    is_unmapped: bool = False,
) -> SimpleNamespace:
    """Build a minimal pysam-AlignedSegment-like object for testing."""
    return SimpleNamespace(
        query_name=query_name,
        reference_start=reference_start,
        query_sequence=query_sequence,
        cigartuples=cigartuples,
        is_reverse=is_reverse,
        is_unmapped=is_unmapped,
    )


# ── build_read_row tests ─────────────────────────────────────────


class TestBuildReadRow:
    """Unit tests for build_read_row CIGAR handling."""

    def test_perfect_match(self):
        """All bases match — every ref position should have aligned base."""
        read = _mock_read("r1", 0, "ACGT", [(CIGAR_MATCH, 4)])
        aligned, inserts = build_read_row(read, 0, 4)
        assert aligned == {0: "A", 1: "C", 2: "G", 3: "T"}
        assert inserts == {}

    def test_seq_match_and_mismatch(self):
        """= (SEQ_MATCH) and X (SEQ_MISMATCH) ops produce aligned bases."""
        read = _mock_read(
            "r1", 0, "ACGT", [(CIGAR_SEQ_MATCH, 2), (CIGAR_SEQ_MISMATCH, 2)]
        )
        aligned, inserts = build_read_row(read, 0, 4)
        assert aligned == {0: "A", 1: "C", 2: "G", 3: "T"}
        assert inserts == {}

    def test_insertion(self):
        """Insertion bases are stored under the anchor position."""
        #  Query:  A C ins(TT) G T
        #  CIGAR:  2M 2I 2M
        read = _mock_read(
            "r1", 0, "ACTTGT", [(CIGAR_MATCH, 2), (CIGAR_INS, 2), (CIGAR_MATCH, 2)]
        )
        aligned, inserts = build_read_row(read, 0, 4)
        assert aligned == {0: "A", 1: "C", 2: "G", 3: "T"}
        assert dict(inserts) == {1: ["T", "T"]}

    def test_deletion(self):
        """Deleted ref positions get '-' in the aligned dict."""
        #  Query:  A C _ _ G T
        #  CIGAR:  2M 2D 2M
        read = _mock_read(
            "r1", 0, "ACGT", [(CIGAR_MATCH, 2), (CIGAR_DEL, 2), (CIGAR_MATCH, 2)]
        )
        aligned, inserts = build_read_row(read, 0, 6)
        assert aligned[0] == "A"
        assert aligned[1] == "C"
        assert aligned[2] == "-"
        assert aligned[3] == "-"
        assert aligned[4] == "G"
        assert aligned[5] == "T"
        assert inserts == {}

    def test_soft_clip(self):
        """Soft-clipped bases are skipped and not placed on the grid."""
        #  Query:  (SS) A C G T (SS)
        #  CIGAR:  2S 4M 2S
        read = _mock_read(
            "r1",
            0,
            "SSACGTSS",
            [(CIGAR_SOFT_CLIP, 2), (CIGAR_MATCH, 4), (CIGAR_SOFT_CLIP, 2)],
        )
        aligned, inserts = build_read_row(read, 0, 4)
        assert aligned == {0: "A", 1: "C", 2: "G", 3: "T"}
        assert inserts == {}

    def test_ref_skip(self):
        """N (REF_SKIP) advances ref without consuming query — gap in aligned."""
        #  Query:  A C ... G T
        #  CIGAR:  2M 3N 2M
        read = _mock_read(
            "r1", 0, "ACGT", [(CIGAR_MATCH, 2), (CIGAR_REF_SKIP, 3), (CIGAR_MATCH, 2)]
        )
        aligned, inserts = build_read_row(read, 0, 7)
        assert aligned[0] == "A"
        assert aligned[1] == "C"
        assert 2 not in aligned
        assert 3 not in aligned
        assert 4 not in aligned
        assert aligned[5] == "G"
        assert aligned[6] == "T"

    def test_window_clipping(self):
        """Only bases within the ref_start..ref_end window are included."""
        read = _mock_read("r1", 0, "ACGTACGT", [(CIGAR_MATCH, 8)])
        aligned, inserts = build_read_row(read, 2, 6)
        assert aligned == {2: "G", 3: "T", 4: "A", 5: "C"}

    def test_read_starting_after_window(self):
        """Read fully outside the window produces empty results."""
        read = _mock_read("r1", 10, "ACGT", [(CIGAR_MATCH, 4)])
        aligned, inserts = build_read_row(read, 0, 5)
        assert aligned == {}
        assert inserts == {}

    def test_insertion_outside_window(self):
        """Insertion anchored outside the window is not recorded."""
        #  Query:  ins(TT) A C G T
        #  CIGAR:  2I 4M at ref pos 0
        #  Insertion anchor = -1 (before window start=0), so ignored.
        read = _mock_read("r1", 0, "TTACGT", [(CIGAR_INS, 2), (CIGAR_MATCH, 4)])
        # Actually this insertion anchor is rpos-1 = 0-1 = -1
        aligned, inserts = build_read_row(read, 0, 4)
        assert aligned == {0: "A", 1: "C", 2: "G", 3: "T"}
        assert inserts == {}

    def test_multiple_insertions_at_same_position(self):
        """Longer insertion produces multiple bases at the anchor."""
        #  CIGAR: 1M 4I 1M at ref 0
        read = _mock_read(
            "r1", 0, "ATTTTC", [(CIGAR_MATCH, 1), (CIGAR_INS, 4), (CIGAR_MATCH, 1)]
        )
        aligned, inserts = build_read_row(read, 0, 2)
        assert aligned == {0: "A", 1: "C"}
        assert dict(inserts) == {0: ["T", "T", "T", "T"]}

    def test_complex_cigar(self):
        """Mixed CIGAR: 3M 2I 1D 2M mimics a realistic read."""
        #  Query:  A C G ins(TT) _ A C
        #  CIGAR:  3M 2I 1D 2M, starting at ref 0
        read = _mock_read(
            "r1",
            0,
            "ACGTTAC",
            [(CIGAR_MATCH, 3), (CIGAR_INS, 2), (CIGAR_DEL, 1), (CIGAR_MATCH, 2)],
        )
        aligned, inserts = build_read_row(read, 0, 6)
        assert aligned[0] == "A"
        assert aligned[1] == "C"
        assert aligned[2] == "G"
        assert aligned[3] == "-"  # deletion
        assert aligned[4] == "A"
        assert aligned[5] == "C"
        assert dict(inserts) == {2: ["T", "T"]}

    def test_bases_uppercased(self):
        """Lowercase query bases are converted to uppercase."""
        read = _mock_read("r1", 0, "acgt", [(CIGAR_MATCH, 4)])
        aligned, inserts = build_read_row(read, 0, 4)
        assert aligned == {0: "A", 1: "C", 2: "G", 3: "T"}


# ── bam_panel integration tests ──────────────────────────────────
# These require pysam and indexed BAM files from examples/.


class TestBamPanel:
    """Integration tests for bam_panel using example files."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_pysam(self):
        pytest.importorskip("pysam")

    def test_indel_bam_panel_shape(self):
        """Panel from indel BAM has correct structure."""
        from tview.bam import bam_panel

        panel = bam_panel("examples/indel_sorted.bam", "examples/ref.fa", "chr1:0-50")
        assert panel.label == "indel_sorted"
        assert len(panel.ref_row) > 0
        assert len(panel.seq_rows) == 5
        assert panel.total_cols >= 50

    def test_indel_bam_insertion_columns(self):
        """Insertion columns are detected in the indel BAM."""
        from tview.bam import bam_panel

        panel = bam_panel("examples/indel_sorted.bam", "examples/ref.fa", "chr1:0-50")
        assert len(panel.ins_columns) > 0

    def test_reads_sorted_by_start(self):
        """Reads in seq_rows are sorted by start position then strand."""
        from tview.bam import bam_panel

        panel = bam_panel("examples/indel_sorted.bam", "examples/ref.fa", "chr1:0-50")
        # read5 is reverse-strand and should come last among reads starting at same position
        names = [name for name, _, _ in panel.seq_rows]
        assert "read5" in names

    def test_col_labels_start_at_one(self):
        """Column labels are 1-based and start at position 1."""
        from tview.bam import bam_panel

        panel = bam_panel("examples/indel_sorted.bam", "examples/ref.fa", "chr1:0-50")
        labels = [lbl for _, lbl in panel.col_labels]
        assert labels[0] == "1"

    def test_sample2_bam(self):
        """Second example BAM loads successfully."""
        from tview.bam import bam_panel

        panel = bam_panel("examples/sample2_sorted.bam", "examples/ref.fa", "chr1:0-50")
        assert panel.label == "sample2_sorted"
        assert len(panel.ref_row) > 0
