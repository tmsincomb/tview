# tview Repository Audit Report

**Date:** 2026-02-28
**Auditor:** Claude (automated)
**Scope:** Full repository audit — code quality, security, CI/CD, testing, documentation

---

## Executive Summary

tview is a well-structured, well-tested Python package for publication-quality alignment visualization. The codebase is clean, modular, and follows Python best practices. No critical security vulnerabilities were found. The findings below are primarily minor improvements and polish items.

**Overall Grade: A-**

| Category | Grade | Notes |
|----------|-------|-------|
| Security | A+ | No vulnerabilities found |
| Code Quality | A | Clean, modular, well-typed |
| Testing | A- | Strong coverage; minor gaps |
| CI/CD | B+ | Functional; minor improvements possible |
| Documentation | B+ | Good README; API docs inaccurate in one spot |

---

## 1. Security Audit

### 1.1 No Vulnerabilities Found

- No use of `eval()`, `exec()`, `pickle`, or unsafe deserialization
- YAML loaded via `yaml.safe_load()` (not `yaml.load()`)
- No subprocess calls, `os.system()`, or shell execution
- No hardcoded credentials or secrets
- File operations use `pathlib.Path` consistently
- CLI input parsed via `click` (safe argument handling)
- `.gitignore` properly excludes `.env`, `.pypirc`, and credential patterns

### 1.2 Dependency Security

Dependencies use minimum version pins (`>=`) rather than locked versions, which is the correct approach for a library — it allows downstream consumers to receive security patches. No known-vulnerable package versions are required.

### 1.3 CI/CD Security

- `publish.yml` uses PyPA trusted publishing (OIDC) — the most secure approach, no API tokens stored
- Minimal permissions model (`contents: read`, `id-token: write`)
- GitHub Actions pinned to major versions (v4, v5)

### 1.4 Low-Risk Observations

- **Integer parsing without bounds checking** (`bam.py:93`, `cli.py:43`): `int()` conversions on user-supplied region/column strings have no upper-bound validation. Risk is negligible since pysam performs its own bounds checking and Python handles arbitrary-precision integers, but defense-in-depth validation could be added.

---

## 2. Code Quality Audit

### 2.1 Strengths

- **Clean modular architecture**: Parsing (bam.py, fasta.py), data model (models.py), rendering (renderer.py), config (config.py), CLI (cli.py) are well-separated
- **Type annotations** throughout via `from __future__ import annotations`
- **Dataclass usage** for `Panel` model — clean, immutable-style design
- **LRU caching** on config loaders prevents redundant file I/O
- **Lazy import** of `pysam` in `bam_panel()` — good for Windows compatibility
- **Consistent formatting** via black (88-char lines) and isort

### 2.2 Issues Found

#### BUG: `GAP_COLOR` imported but never used in renderer.py

`config.py` exports `GAP_COLOR` and it's imported in `test_config.py`, but `renderer.py:14-27` imports many config constants yet `GAP_COLOR` is not among them. In `renderer.py:191` and `renderer.py:215`, gap characters (`"-"`) are rendered using `TEXT_COLOR` (#000000) rather than `GAP_COLOR` (#9E9E9E). This means gaps render as black dashes rather than the intended grey color.

**File:** `renderer.py:191,215`
**Severity:** Low (visual bug)
**Fix:** Import `GAP_COLOR` and use it for gap rendering instead of `TEXT_COLOR`.

#### BUG: `mono` and `mono_sm` are identical in `_resolve_font()`

`renderer.py:50-52` creates two `FontProperties` objects (`mono` and `mono_sm`) that are always identical — same font path, same size, same weight. The naming suggests `mono_sm` should be a smaller font for tick labels, but it's the same size as `mono`.

**File:** `renderer.py:50-52, 63-64`
**Severity:** Low (cosmetic — tick labels should arguably be smaller)
**Fix:** Either make `mono_sm` use a smaller `fontsize` (e.g., `fontsize * 0.8`) or remove the distinction and use a single `FontProperties` object.

#### ISSUE: `render_panels()` prints to stdout unconditionally

`renderer.py:328` calls `print(f"Saved: {out_path} ...")` on every render. This is fine for CLI usage but problematic when used as a library — callers cannot suppress this output without redirecting stdout.

**File:** `renderer.py:328`
**Severity:** Low
**Fix:** Use `logging.info()` instead, or remove the print and let the CLI handle status messages.

#### ISSUE: Mutable default in `Panel.ins_columns`

`models.py:41` uses `None` as default for `ins_columns` and converts to `set()` in `__post_init__`. This is correctly done (avoids the mutable default pitfall), but the type annotation `set[int] | None = None` means callers see `Optional[set]` in their IDE. Since `__post_init__` always normalizes to `set()`, the external type should be `set[int]` with `field(default_factory=set)`.

**File:** `models.py:41-45`
**Severity:** Very low (type hint accuracy)

#### ISSUE: `fasta_panel` doesn't validate equal sequence lengths

`fasta.py:92-132`: When reading aligned FASTA files, there's no check that all sequences have equal length. If sequences differ in length, `row += ["-"] * (aln_len - len(row))` silently pads shorter sequences, but sequences longer than the reference are silently truncated via `list(seq.upper()[:aln_len])`. This could mask input errors.

**File:** `fasta.py:111`
**Severity:** Low
**Fix:** Add an optional warning when sequence lengths differ from the reference.

#### ISSUE: Homepage URL mismatch

`pyproject.toml:45-47` lists `Homepage` and `Repository` as `https://github.com/MurrellGroup/tview`, but the repo remote and image URLs in `README.md:7-17` point to `https://github.com/tmsincomb/tview`. These should be consistent.

**File:** `pyproject.toml:45-47`, `README.md:7-17`
**Severity:** Low (broken links after any transfer)

---

## 3. Testing Audit

### 3.1 Strengths

- **~60 test cases** across 3 test modules with good coverage
- **Property-based testing** via Hypothesis — excellent for edge cases
- **Visual regression artifacts** saved to `tests/output/` for manual inspection
- **Cross-platform CI** (Linux, macOS, Windows) x (Python 3.9, 3.13)
- **patchworklib integration tests** with proper `skipif` guard
- **CLI tested via `click.testing.CliRunner`** — best practice

### 3.2 Gaps

#### MISSING: No dedicated `test_bam.py`

BAM parsing (`bam.py`) is only tested indirectly via `test_cli.py:test_mixed_bam_and_fasta`. There are no unit tests for:
- `build_read_row()` with various CIGAR operations
- Insertion column expansion logic
- Edge cases: reads entirely outside the window, zero-length regions, overlapping reads

**Severity:** Medium — BAM CIGAR parsing is the most complex logic in the codebase and deserves dedicated unit tests.

#### MISSING: No test for `read_fasta` with non-existent file

The docstring for `read_fasta` shows a `FileNotFoundError` doctest, but there's no pytest test for this error path.

#### MISSING: No test for `fasta_panel` with empty file

`fasta_panel` raises `ValueError` for empty files, but there's no test exercising this path.

#### MISSING: No test for `_resolve_font()` fallback paths

Font resolution has multiple fallback branches (`renderer.py:45-65`) but no tests for when preferred fonts aren't available.

#### MINOR: Test output directory committed to git?

`conftest.py` creates `tests/output/` and tests write artifacts there. This directory doesn't appear in `.gitignore`. If test artifacts are committed, they bloat the repo; if not, the `mkdir` in `conftest.py` handles it at runtime.

**Fix:** Add `tests/output/` to `.gitignore` to prevent accidental commits of rendered test images.

---

## 4. CI/CD Audit

### 4.1 Strengths

- **Cross-platform matrix**: 3 OS x 2 Python versions = 6 job combinations
- **`fail-fast: false`**: All matrix jobs run even if one fails — good for debugging
- **Trusted PyPI publishing**: Modern OIDC-based approach, no stored secrets
- **`setuptools-scm`**: Version derived from git tags, no manual version bumps

### 4.2 Issues

#### ISSUE: CI tests only Python 3.9 and 3.13, but classifiers list 3.9-3.14

`pyproject.toml` classifies support for Python 3.9, 3.10, 3.11, 3.12, 3.13, and 3.14, but CI only tests 3.9 and 3.13. This is a reasonable compromise for cost, but consider adding at least 3.12 (current stable) to catch version-specific issues.

**File:** `.github/workflows/ci.yml:15`
**Severity:** Low

#### ISSUE: `publish.yml` uses `if: always() && needs.build.result == 'success'`

`publish.yml:46`: The `always()` condition with a success check is unusual. A simpler `if: success()` (or just the default behavior) would suffice since `needs: [build]` already gates on the build job.

**File:** `.github/workflows/publish.yml:46`
**Severity:** Very low

#### ISSUE: No linting in CI

CI runs `pytest` but doesn't run `black --check` or `isort --check`. Pre-commit hooks exist locally but aren't enforced in CI, meaning PRs can merge with formatting violations.

**File:** `.github/workflows/ci.yml`
**Severity:** Low
**Fix:** Add a lint step: `black --check src/ tests/` and `isort --check src/ tests/`.

#### ISSUE: GitHub Actions not pinned to commit SHAs

Actions are pinned to major versions (`actions/checkout@v4`) rather than full commit SHAs. While this is common practice, SHA pinning (e.g., `actions/checkout@<sha>`) is more secure against supply chain attacks.

**Severity:** Very low (standard practice for most projects)

---

## 5. Documentation Audit

### 5.1 Strengths

- **Comprehensive README** with visual examples, CLI reference, Python API, and publication tips
- **Docstrings with doctests** on all public functions
- **Visual convention table** clearly explains rendering symbols
- **Color palette documentation** with hex values
- **Stdin piping documentation** — good UX touch

### 5.2 Issues

#### BUG: Python API example uses wrong parameter names

`README.md:150` shows:
```python
panel = fasta_panel("aligned.fasta", col_start=1, col_end=120)
```

But `fasta_panel()` actually accepts `columns: list[int] | None`, not `col_start`/`col_end`. The correct call would be:
```python
panel = fasta_panel("aligned.fasta", columns=list(range(1, 121)))
```

Same issue at `README.md:188`.

**Severity:** Medium (users copying this example will get a `TypeError`)

#### ISSUE: `--columns` description inconsistency

README line 272 says `--columns TEXT  Column range for FASTA, 1-based inclusive (e.g. 1-120)` in the argument reference section, but the CLI help (line 92) says `Column positions for FASTA, 1-based (e.g. 1-120, 5,40,690, or 5,10-20,40)`. The CLI help is more accurate since `--columns` now supports discrete positions and mixed ranges (added in recent commit), but the README argument table hasn't been fully updated.

**Severity:** Low

#### MISSING: No `CONTRIBUTING.md` or development setup guide

For an open-source project, adding basic contribution guidelines (how to run tests, formatting requirements, PR process) would help new contributors.

**Severity:** Very low

#### ISSUE: `patchworklib` install extra doesn't exist

`README.md:180` shows `pip install tview[compose]` but `pyproject.toml` doesn't define a `[compose]` optional dependency group. `patchworklib` is already a core dependency in `dependencies`, so this extra is unnecessary, but the documentation is misleading.

**Severity:** Low (confusing but non-breaking since patchworklib is already a dependency)

---

## 6. Architecture Observations

### 6.1 Good Design Decisions

- **Separation of `draw_panels()` and `render_panels()`**: Enables both quick CLI usage and advanced composition with patchworklib/matplotlib
- **YAML-driven configuration**: Colors and fonts are externalizable without code changes
- **`Panel` as the interchange format**: Clean data boundary between parsers and renderer
- **Lazy pysam import**: Enables the package to install and run FASTA-only workflows on Windows

### 6.2 Potential Improvements

- **Consider using `logging` module**: Replace `print()` in `renderer.py:328` and potential future debug output with proper logging levels
- **Consider adding a `U` base to nucleotide palette**: For RNA sequence support (`U` maps to the same role as `T`)
- **Consider adding `--width`/`--height` CLI options**: To override automatic figure sizing for specific publication requirements

---

## 7. Summary of Actionable Findings

### Must Fix (Bugs)

| # | File | Issue | Severity |
|---|------|-------|----------|
| 1 | `README.md:150,188` | Python API examples use `col_start`/`col_end` — parameters don't exist | Medium |
| 2 | `renderer.py:191,215` | Gaps use `TEXT_COLOR` (black) instead of `GAP_COLOR` (grey) | Low |

### Should Fix

| # | File | Issue | Severity |
|---|------|-------|----------|
| 3 | Tests | No dedicated `test_bam.py` for CIGAR parsing logic | Medium |
| 4 | `ci.yml` | No linting step (`black --check`, `isort --check`) | Low |
| 5 | `pyproject.toml:45-47` | Homepage URL points to `MurrellGroup/tview`, not `tmsincomb/tview` | Low |
| 6 | `README.md:180` | `pip install tview[compose]` extra doesn't exist | Low |
| 7 | `renderer.py:328` | `print()` should be `logging.info()` for library use | Low |

### Nice to Have

| # | File | Issue | Severity |
|---|------|-------|----------|
| 8 | `renderer.py:50-52` | `mono_sm` is identical to `mono` — should differ | Very Low |
| 9 | `.gitignore` | Add `tests/output/` to prevent accidental commits | Very Low |
| 10 | `ci.yml:15` | Add Python 3.12 to CI matrix | Very Low |
| 11 | `models.py:41` | Use `field(default_factory=set)` instead of `None` + `__post_init__` | Very Low |
