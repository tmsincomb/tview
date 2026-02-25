"""Shared fixtures for tview tests."""

from __future__ import annotations

from pathlib import Path

import pytest

OUTPUT_DIR = Path(__file__).parent / "output"


@pytest.fixture(scope="session", autouse=True)
def ensure_output_dir():
    """Create tests/output/ once per session."""
    OUTPUT_DIR.mkdir(exist_ok=True)


@pytest.fixture
def output_dir():
    """Persistent output directory for visual inspection."""
    return OUTPUT_DIR


@pytest.fixture
def write_fasta(tmp_path):
    """Factory fixture that writes sequences to a FASTA file and returns the path."""

    def _write(sequences: list[tuple[str, str]], filename: str = "test.fasta") -> Path:
        path = tmp_path / filename
        with open(path, "w") as f:
            for name, seq in sequences:
                f.write(f">{name}\n{seq}\n")
        return path

    return _write
