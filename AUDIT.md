# tview Package Audit

**Date:** 2025-02-25
**Scope:** Full audit of the `tview` package — code quality, correctness, security, API design, concurrency, error handling, tests, performance, and packaging.

---

## Executive Summary

`tview` is a well-structured, focused Python package (~650 lines across 3 source files) for rendering publication-quality alignment visualizations from BAM and FASTA files. The code is clean, the test suite uses property-based testing (Hypothesis), and the CLI is ergonomic. However, there are several issues ranging from bugs to design concerns that should be addressed before a stable 1.0 release.

**Severity legend:** 🔴 High — 🟡 Medium — 🟢 Low — ℹ️ Informational

---

## 1. Correctness & Bugs

### 🔴 B01: Ambiguous region coordinate convention (`tview.py:301-302`)

The `bam_panel` region parser passes user coordinates directly to pysam, which uses **0-based half-open** intervals. However, the bioinformatics convention for region strings like `chr1:100-200` is **1-based inclusive** (as used by samtools, IGV, UCSC, etc.). The README examples show `chr1:100-200` without clarifying the coordinate system, which will confuse users.

```python
# Current: 0-based half-open (pysam convention)
start, end = [int(x) for x in rest.replace(",", "").split("-")]
```

**Recommendation:** Either (a) convert from 1-based inclusive to 0-based half-open internally (subtract 1 from start), or (b) clearly document that coordinates are 0-based. Option (a) matches user expectations.

### 🟡 B02: `mono_sm` is identical to `mono` (`tview.py:391-393, 407-408`)

`_resolve_font` creates `mono_sm` with the exact same size and weight as `mono`. The variable name suggests it should be smaller (for tick labels). This is dead differentiation — either the intended size difference was lost, or the variable should be consolidated.

```python
mono = fm.FontProperties(fname=helv_path, size=fontsize, weight="bold")
mono_sm = fm.FontProperties(fname=helv_path, size=fontsize, weight="bold")  # identical
```

### 🟡 B03: `render_panels` uses pyplot global state — not thread-safe (`tview.py:451, 566-575`)

The function uses `plt.subplots()`, `plt.subplots_adjust()`, `plt.savefig()`, and `plt.close()`. If an exception occurs mid-render, the figure leaks. Even without exceptions, the pyplot state machine is not thread-safe.

```python
# Current (global state)
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
plt.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.01)
plt.savefig(out_path, ...)
plt.close()
```

**Recommendation:** Use the OO interface with a `try/finally` for cleanup:
```python
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
try:
    ...
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.01)
    fig.savefig(out_path, ...)
finally:
    plt.close(fig)
```

### 🟡 B04: `fasta_panel` column labels reset to 1 after slicing (`tview.py:222-229`)

When `col_start` and `col_end` are used, the reference position counter restarts at 1 for the sliced window. If a user requests columns 50–100, the x-axis shows "1, 10, 20..." instead of positions relative to the full alignment or the original reference numbering. This could be confusing.

### 🟡 B05: Missing CIGAR hard clip constant (`tview.py:82-89`)

CIGAR operation 5 (hard clip, `H`) and 6 (padding, `P`) are not defined. They are safely ignored by the if/elif chain (they consume neither query nor reference bases), but defining them as constants with a comment would improve clarity and prevent confusion during maintenance.

### 🟢 B06: Potential `None` query_sequence in `build_read_row` (`tview.py:264`)

`read.query_sequence` can be `None` for some pysam reads (e.g., reads with no stored sequence). Accessing `read.query_sequence[qpos]` would raise a `TypeError`. The `bam_panel` filter (`r.cigartuples` being truthy) correlates with having a sequence, but doesn't guarantee it.

### 🟢 B07: Outdated module docstring (`tview.py:3`)

The module docstring says `tview_image.py` — this is a leftover from before the package was restructured.

---

## 2. Security

### 🟢 S01: No input path validation for `--bam` and `--fasta`

`--ref` uses `click.Path(exists=True)` for validation, but `--bam` and `--fasta` accept arbitrary strings without existence checks. This means users get a raw Python traceback (FileNotFoundError or pysam error) rather than a helpful Click error message.

### ℹ️ S02: Security surface is minimal

As a local CLI tool, the attack surface is negligible. The tool reads files the user specifies and writes images. There are no network operations, no deserialization of untrusted data, and no shell command construction. If used as a library in a web service context, file path inputs should be sanitized by the caller.

---

## 3. API Design

### 🔴 A01: `matplotlib.use("Agg")` at import time (`tview.py:37`)

Calling `matplotlib.use("Agg")` at module level permanently changes the matplotlib backend for the entire process. This breaks interactive backends (Qt, Tk, etc.) for any library consumer who imports `tview` alongside interactive matplotlib usage.

**Recommendation:** Move backend selection into `render_panels` or guard it:
```python
import matplotlib
if matplotlib.get_backend() == "":
    matplotlib.use("Agg")
```

### 🟡 A02: `render_panels` prints to stdout (`tview.py:576`)

```python
print(f"Saved: {out_path} ({dpi} dpi, {len(panels)} panel(s), " f"{max_cols} cols)")
```

A library function should not print to stdout. CLI consumers can print their own messages; library consumers may want silent operation.

**Recommendation:** Use `logging.info()` or add a `quiet` parameter, or move the print to the CLI layer.

### 🟡 A03: Inconsistent coordinate conventions between FASTA and BAM modes

- `fasta_panel`: 1-based inclusive `col_start`/`col_end`
- `bam_panel`: 0-based half-open via raw passthrough to pysam

Users switching between modes will be confused by this inconsistency.

### 🟢 A04: `Panel.ins_columns` default handling (`tview.py:120-124`)

Using `set[int] | None = None` with `__post_init__` conversion is fine but could be simplified with `field(default_factory=set)` to avoid the `None` intermediate state entirely.

### 🟢 A05: `render_panels` returns `None`

The function could return the output path or the figure object for further programmatic use (e.g., compositing).

---

## 4. Error Handling & Edge Cases

### 🟡 E01: No validation for empty `panels` list (`tview.py:440`)

```python
max_cols = max(p.total_cols for p in panels)  # ValueError if panels is empty
```

**Recommendation:** Add `if not panels: raise ValueError("At least one panel is required")`.

### 🟡 E02: Malformed region string produces unhelpful errors (`tview.py:301-302`)

```python
chrom, rest = region.split(":")  # ValueError if no ":"
start, end = [int(x) for x in rest.replace(",", "").split("-")]  # ValueError if no "-"
```

**Recommendation:** Wrap with a clear error message:
```python
try:
    chrom, rest = region.split(":")
    start, end = [int(x) for x in rest.replace(",", "").split("-")]
except (ValueError, IndexError):
    raise ValueError(f"Invalid region format '{region}'. Expected 'chrom:start-end'.")
```

### 🟡 E03: Malformed `--columns` produces unhelpful errors (`cli.py:91-93`)

Same issue as E02 — invalid input like `--columns abc` raises an unhandled `ValueError`.

### 🟢 E04: Single-sequence FASTA produces a panel with 0 seq_rows

`fasta_panel` accepts a file with only 1 sequence (the reference) and creates a panel with `seq_rows=[]`. This renders as just a reference row, which may confuse users expecting reads.

### 🟢 E05: FASTA sequences longer than reference are silently truncated (`tview.py:218`)

```python
row = list(seq.upper()[:aln_len])  # truncates to ref length
```

This is reasonable for aligned FASTA (all seqs should be the same length) but could mask alignment file corruption.

---

## 5. Performance

### 🟡 P01: Individual text artists per base (`tview.py:474-533`)

Each base position creates a separate `ax.text()` call and potentially an `ax.add_patch()` for mismatches. For a 200-column × 50-sequence alignment, this creates ~10,000 text artists and potentially thousands of patches. Rendering becomes very slow for large alignments.

**Recommendation:** For larger alignments, consider using `PatchCollection` for batch rectangle rendering and investigate whether a rasterized approach (pixel grid) might be more appropriate above a certain size threshold.

### ℹ️ P02: No alignment size limits

There's no guardrail preventing a user from attempting to render a 10,000-column × 1,000-sequence alignment, which would exhaust memory and time. A warning for large inputs would be helpful.

---

## 6. Test Quality

### 🟡 T01: No unit tests for `build_read_row`

The core CIGAR parsing logic in `build_read_row` is only tested indirectly through the `test_mixed_bam_and_fasta` integration test. Given its complexity (handling M, I, D, N, S, =, X operations), it deserves dedicated unit tests with known CIGAR strings.

### 🟡 T02: Limited BAM test coverage

Only one BAM test (`test_mixed_bam_and_fasta`) exists, and it uses `importorskip`. There are no tests for:
- Reads with insertions/deletions
- Reverse-strand reads
- Reads that partially overlap the region
- Edge cases like empty BAM regions

### 🟢 T03: Temp file leak risk in `_write_fasta` (`test_fasta.py:20-26`)

Uses `NamedTemporaryFile(delete=False)`. While tests clean up in `finally` blocks, the helper itself doesn't ensure cleanup if used carelessly. The `conftest.py` `write_fasta` fixture is a better pattern.

### 🟢 T04: Class-level mutable counter in tests (`test_fasta.py:168-169`)

`TestRenderPanels._nt_counter` and `_aa_counter` are class-level state modified in Hypothesis tests. This could produce non-deterministic file names with parallel test execution (pytest-xdist).

### ℹ️ T05: Good use of Hypothesis

The property-based testing with `aligned_nt_seqs` and `aligned_aa_seqs` strategies is well done and catches edge cases that manual tests would miss.

---

## 7. Packaging & CI

### 🟡 PKG01: `_version.py` is tracked in git

The file header says "don't track in version control," but it's committed. This causes unnecessary merge conflicts and dirty working trees when building locally. It should be added to `.gitignore`.

### 🟢 PKG02: No linting or type checking in CI (`.github/workflows/ci.yml`)

Pre-commit hooks exist for `black` and `isort`, but CI doesn't run them. There's also no `mypy` or `ruff` for static analysis. CI should enforce the same standards as pre-commit.

### 🟢 PKG03: Python 3.14 classifier may be premature

`Programming Language :: Python :: 3.14` is listed in classifiers. Python 3.14 may not be fully released or tested against.

### 🟢 PKG04: Publish workflow `if` condition (`publish.yml:45`)

```yaml
if: always() && needs.build.result == 'success'
```

Using `always()` means this job runs even if the workflow is cancelled. A safer pattern is `if: needs.build.result == 'success'` (without `always()`), which naturally skips on cancellation.

---

## 8. Documentation

### 🟢 D01: Region format not documented

The `--region` help text says `chr:start-end` but doesn't specify 0-based vs 1-based. The README examples also don't clarify. This is critical for bioinformatics users who expect 1-based.

### 🟢 D02: `--columns` with single value undocumented

`--columns 50` (without end) is parsed as `col_start=50, col_end=None`, meaning "columns 50 to end." This behavior is not documented.

### ℹ️ D03: Good README overall

The README is comprehensive with clear examples, visual convention tables, and practical tips for publication figures.

---

## Summary of Recommendations by Priority

| Priority | Count | Key Items |
|----------|-------|-----------|
| 🔴 High | 2 | Region coordinate convention (B01), matplotlib backend at import (A01) |
| 🟡 Medium | 12 | Thread safety (B03), stdout printing (A02), empty panels crash (E01), missing CIGAR tests (T01), _version.py in git (PKG01) |
| 🟢 Low | 11 | Path validation (S01), font duplication (B02), minor edge cases |
| ℹ️ Info | 4 | Good test practices, low security surface, good README |

---

## Files Reviewed

| File | Lines | Description |
|------|-------|-------------|
| `src/tview/tview.py` | 657 | Core library — panel builders, renderer |
| `src/tview/cli.py` | 99 | Click CLI entry point |
| `src/tview/__init__.py` | 27 | Public API exports |
| `src/tview/_version.py` | 34 | Auto-generated version (setuptools-scm) |
| `tests/test_fasta.py` | 289 | FASTA and rendering tests (Hypothesis) |
| `tests/test_cli.py` | 132 | CLI integration tests |
| `tests/conftest.py` | 36 | Test fixtures |
| `pyproject.toml` | 67 | Build configuration |
| `.github/workflows/ci.yml` | 37 | CI workflow |
| `.github/workflows/publish.yml` | 61 | PyPI publish workflow |
