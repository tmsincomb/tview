# TODO

## patchworklib integration

**Status:** Suggested but not implemented or tested.

- Listed as dependency in `pyproject.toml:39` (`patchworklib>=0.6`).
- Not imported anywhere in `src/tview/`.
- Mentioned in `src/tview/renderer.py:94,149` docstrings only — `draw_panels` claims compatibility with `patchworklib.Brick` but no code path exercises it.
- Tests exist at `tests/test_fasta.py:558-625` (`TestPatchworklib`) but skip whenever `patchworklib` is not installed in the active env. Currently skipped under the `tview` conda env.

### Action items

- [ ] Install `patchworklib` in the `tview` conda env so `TestPatchworklib` runs.
- [ ] Verify `draw_panels` actually renders onto `pw.Brick` objects (claimed in docstring).
- [ ] Add a real example in `examples/` or `README.md` demonstrating multi-panel composition with patchworklib.
- [ ] Decide whether `patchworklib` should remain a hard dependency or move to an optional extras group (e.g. `[project.optional-dependencies] compose = ["patchworklib>=0.6"]`).
