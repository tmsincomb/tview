"""Tests for the Click CLI."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from tview.cli import main

from .conftest import OUTPUT_DIR

EXAMPLES = Path(__file__).parent.parent / "examples"


class TestCLI:
    def test_help(self):
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "Publication-quality" in result.output

    def test_no_args_errors(self):
        runner = CliRunner()
        result = runner.invoke(main, [])
        assert result.exit_code != 0

    def test_fasta_mode(self, write_fasta):
        seqs = [("ref", "ACGTACGT"), ("s1", "ACCTACGT")]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_fasta_mode.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0
        assert out.exists()

    def test_fasta_with_columns(self, write_fasta):
        seqs = [("ref", "ACGTACGTACGT"), ("s1", "ACCTACGTACGT")]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_columns_1-8.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "--columns",
                "1-8",
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0
        assert out.exists()

    def test_fasta_with_discrete_columns(self, write_fasta):
        seqs = [("ref", "ACGTACGTACGT"), ("s1", "ACCTACGTACGT")]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_columns_discrete.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "--columns",
                "1,3,5,7",
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0
        assert out.exists()

    def test_fasta_with_mixed_columns(self, write_fasta):
        seqs = [("ref", "ACGTACGTACGT"), ("s1", "ACCTACGTACGT")]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_columns_mixed.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "--columns",
                "1-4,8,12",
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0
        assert out.exists()

    def test_bam_without_ref_errors(self):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--bam",
                "fake.bam",
                "--region",
                "chr1:1-10",
            ],
        )
        assert result.exit_code != 0

    def test_mixed_bam_and_fasta(self, write_fasta):
        """BAM + FASTA panels stacked in one figure."""
        pysam = pytest.importorskip("pysam")  # noqa: F841
        fasta_seqs = [
            ("chr1_ref", "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("sample_x", "CGATCAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
        ]
        fasta_path = write_fasta(fasta_seqs, "mixed_fasta.fasta")
        out = OUTPUT_DIR / "mixed_bam_fasta.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--bam",
                str(EXAMPLES / "indel_sorted.bam"),
                "--ref",
                str(EXAMPLES / "ref.fa"),
                "--region",
                "chr1:0-50",
                "--fasta",
                str(fasta_path),
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_show_invokes_kitten(self, write_fasta, monkeypatch):
        calls = []
        monkeypatch.setattr(
            "tview.cli.shutil.which",
            lambda cmd: "/fake/bin/kitten" if cmd == "kitten" else None,
        )
        monkeypatch.setattr(
            "tview.cli.subprocess.run",
            lambda cmd, **kw: calls.append(cmd) or None,
        )
        seqs = [("ref", "ACGT"), ("s1", "ACCT")]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_show_invokes.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "--show",
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert len(calls) == 1
        assert calls[0][0] == "/fake/bin/kitten"
        assert "icat" in calls[0]
        assert str(out) in calls[0]

    def test_show_missing_kitten_is_soft_failure(self, write_fasta, monkeypatch):
        monkeypatch.setattr("tview.cli.shutil.which", lambda cmd: None)
        seqs = [("ref", "ACGT"), ("s1", "ACCT")]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_show_missing.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "--show",
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0
        assert out.exists()
        assert "kitten" in result.stderr

    def test_max_rows_limits_fasta(self, write_fasta):
        seqs = [("ref", "ACGT")] + [(f"s{i}", "ACCT") for i in range(10)]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_max_rows.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "--max-rows",
                "3",
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_palette_aa(self, write_fasta):
        seqs = [("ref", "MRVKGSTN"), ("s1", "MRVKDEWF")]
        fasta_path = write_fasta(seqs)
        out = OUTPUT_DIR / "cli_palette_aa.png"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--fasta",
                str(fasta_path),
                "--palette",
                "aa",
                "-o",
                str(out),
                "--dpi",
                "150",
            ],
        )
        assert result.exit_code == 0
        assert out.exists()
