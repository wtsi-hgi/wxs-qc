# AGENTS.md

Instructions for AI coding agents working in this repository.

## Project Context

WxS-QC is a quality-control pipeline for human germline short-variant WGS and WES cohorts.
It is based on Hail and draws from gnomAD QC pipelines,
with a number of local methods and workflow conventions layered on top.

This is an old and complicated codebase.
It contains numbered pipeline stages, standalone scripts, shared helper modules, Hail/Spark-specific behavior, docs,
and tests that cover only a limited part of the actual behavior.
Treat the existing code as the source of truth, even when it differs from the current
development guidelines or from what a new implementation would ideally look like.

The current development direction is iterative improvement of coding style,
maintainability, and reliability.
Changes should be small, scoped, and easy to review.

## Repository Structure

- `0-resource_preparation/`: numbered scripts for preparing external resources.
- `1-import_data/`: import and validation scripts for input data and annotations.
- `2-sample_qc/`: sample QC pipeline steps.
- `3-variant_qc/`: variant QC, random forest, ranking, binning, plotting, and
  filtering steps.
- `4-genotype_qc/`: genotype hard-filter evaluation, VCF export, and downstream
  summaries.
- `wes_qc/`: main shared Python package for pipeline utilities and reusable
  logic.
- `utils/`: older or auxiliary utility code still used by parts of the pipeline.
- `config/`: example/public configuration.
- `docs/`: user and developer documentation.
- `scripts/`: shell helpers, notebooks, and cluster/local execution utilities.
- `tests/`: unit and integration tests. Coverage is limited and should not be
  interpreted as proving that broad refactors are safe.

## Non-Negotiable Scope Rules

- Only modify files that are clearly inside the approved task scope.
- Never touch unrelated pipeline stages, helper modules, docs, config, tests, or
  formatting outside the approved scope.
- Do not perform opportunistic cleanup, broad formatting, import sorting, or
  refactoring in files that are not part of the requested change.
- If a required fix appears to need changes outside the approved scope,
  stop and ask the user for confirmation before editing those files.
- If task instructions conflict with the actual codebase, stop and ask for
  further instructions. Do not guess which side should win.
- If the README, docs, tests, type hints, comments, or function names disagree
  with executable behavior, inspect the code and ask before making assumptions that affect behavior.
- Preserve user changes. Do not revert, overwrite, or reformat work that was not part of your own change.

## Environment Assumptions

- Assume the required configuration, tools, credentials, permissions, and data access have been set up correctly for the assigned task.
- NEVER work around missing tooling, unavailable data, permission errors, cloud access issues,
  or broken local configuration by changing project code or
  widening scope.
- If tooling, permissions, credentials, environment variables, data access,
  or local configuration appear to be missing or incorrect, stop and ask the user to fix the environment before continuing.

## Working Style

- Read the relevant scripts and helper functions before editing. This repository
  has historical patterns and implicit contracts that are not always documented.
- Prefer the existing local style in the file being changed over introducing a
  new pattern.
- Keep changes minimal and behavior-preserving unless the task explicitly asks
  for behavior changes.
- When improving style or reliability, do it incrementally around the requested
  code path. Avoid large architectural changes.
- Use structured parsing APIs for config, YAML, Python, or tabular data when
  practical. Avoid fragile ad hoc text manipulation for code or data formats.
- Add comments only where they clarify non-obvious Hail/Spark behavior,
  historical constraints, or complex transformations.
- Do not rename public scripts, pipeline step files, config keys, or output paths
  unless the user explicitly approves the migration.

## Python and Hail Conventions

- The project targets Python 3.12, but some parts of the code may still assume Python 3.9 compatibility.
  Convert them to the new Python 3.12 syntax when they are in scope.
- Hail MatrixTable/Table reads and writes are expensive and often tied to workflow state.
  Do not move IO across function boundaries unless the task explicitly requires it.
- Prefer keeping Hail data loading and saving in `main()`-style script entry points when that matches nearby code.
- Pipeline step functions should generally accept and return Hail objects rather
  than reading and writing them internally.
- Preserve numbered script sequencing and command-line interfaces unless the
  requested task is specifically about changing them.
- Be careful with path handling. Convert paths for Spark/Hail only at the point where Spark/Hail APIs need that format.

## Tests and Verification

Current automated test coverage is limited and not comprehensive.
The unit tests are currently broken and are not used for validation.
Use the integration smoke tests as the active validation path; they mainly check
for the absence of unhandled exceptions.
Passing integration smoke tests does not guarantee that a change is safe across
the whole pipeline.

Available commands include:

```bash
make integration-test-trios
make integration-test-non-trios
make test-it-one-step test=test_trios_<step-name-or-prefix>
pre-commit run --files <changed-files>
pre-commit run mypy --hook-stage manual
```

Notes:

- Run `make integration-test-trios` to test the pipeline end-to-end through the
  trio integration smoke-test path.
- Run `make test-it-one-step test=test_trios_...` to individually test and debug
  a specific trio pipeline step.
- Use `make integration-test-non-trios` only when the approved change
  specifically affects the non-trio path.
- The current Makefile does not define a plain `integration-test` target, and
  `make test` depends on that missing target. Prefer the explicit trio/non-trio
  integration targets until the Makefile is corrected.
- Integration tests may require Hail/Spark setup and generated intermediate files.
- `PYSPARK_PYTHON` and `PYSPARK_DRIVER_PYTHON` are set by the Makefile.
- Prefer running the smallest relevant test target first.
- If tests cannot be run because of environment, data, cloud, Spark, or time
  constraints, state that clearly in the final response.
- When adding or changing behavior, add focused tests if the touched area has an
  appropriate existing test pattern. Do not build a large new test framework as
  part of a narrow task.

## Documentation

- Keep documentation changes scoped to the behavior or workflow being changed.
- Do not rewrite broad README or docs sections unless requested.
- If code behavior and documentation disagree, stop and ask before changing one
  to match the other, unless the task explicitly identifies which source should
  be corrected.

## Final Response Expectations

When finishing a task, summarize:

- What changed.
- Which files were touched.
- Which tests or checks were run.
- Any tests or checks that were not run and why.
- Any out-of-scope issue discovered that still needs user approval.
