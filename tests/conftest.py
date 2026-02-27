"""Shared fixtures for tview tests."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest

OUTPUT_DIR = Path(__file__).parent / "output"


@pytest.fixture(scope="session", autouse=True)
def ensure_output_dir() -> None:
    """Create tests/output/ once per session."""
    OUTPUT_DIR.mkdir(exist_ok=True)


@pytest.fixture
def output_dir() -> Path:
    """Persistent output directory for visual inspection."""
    return OUTPUT_DIR


@pytest.fixture
def write_fasta(tmp_path: Path) -> Callable[..., Path]:
    """Factory fixture that writes sequences to a FASTA file and returns the path."""

    def _write(sequences: list[tuple[str, str]], filename: str = "test.fasta") -> Path:
        path = tmp_path / filename
        lines = [f">{name}\n{seq}\n" for name, seq in sequences]
        path.write_text("".join(lines))
        return path

    return _write
