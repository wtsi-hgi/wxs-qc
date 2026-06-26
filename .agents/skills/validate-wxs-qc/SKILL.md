---
name: validate-wxs-qc
description: Validate WxS-QC changes by running integration smoke tests, reviewing pipeline contracts, checking scope discipline, and documenting blocked integration/Hail/Spark validation.
---

# Validate WxS-QC Skill

Use this skill when asked to review or validate implementation quality.

## Workflow

1. Run `.agents/skills/validate-wxs-qc/scripts/validate.sh` for baseline checks when practical.
2. Read `AGENTS.md`.
3. Read `artifacts/1_plan.md` and `artifacts/2_implement.md` when present.
4. Review changed files against the approved scope and actual codebase behavior.
5. Check Hail/Spark behavior, config compatibility, numbered script sequence, CLI arguments, output paths, and shared helper call sites affected by the change.
6. Save the validation report to `artifacts/3_validate.md` when using the artifacts workflow.
7. Report findings by severity with file/line references.
8. If no findings, state that explicitly and note residual risks.

## Required Checks

Use the integration smoke tests as the validation baseline:

- `pre-commit run --files <changed-files>`
- `pre-commit run mypy --hook-stage manual` when type changes justify it
- `make test-it-one-step test=test_trios_...` to individually test and debug a specific trio pipeline step
- `make integration-test-trios` to test the pipeline end-to-end through the trio smoke-test path
- `make integration-test-non-trios` only when the approved change specifically affects the non-trio path

Notes:

- Unit tests are currently broken and not used for validation.
- Integration tests may require Hail/Spark setup and generated intermediate files.
- The current Makefile does not define a plain `integration-test` target, and `make test` depends on that missing target.
- Treat blocked checks as residual risk unless equivalent evidence exists.
- Passing integration smoke tests do not prove broad pipeline safety because coverage is limited.

## Review Focus

- Scope matches the approved plan and `AGENTS.md`.
- No unrelated files, stages, docs, config, or formatting were changed.
- Public scripts, CLI flags, config keys, and output paths remain compatible unless explicitly approved.
- Hail Table/MatrixTable schemas, keys, annotations, and IO behavior are preserved or intentionally changed.
- Spark path conversion happens only where needed.
- Shared helper changes have known call sites reviewed.
- Integration smoke tests or manual validation cover the changed behavior where practical.

## Reporting Format

Start with:

Status: PASS / FAIL / BLOCKED

Then list findings ordered by severity:

High:
- Likely bug, regression, data loss, invalid QC result, or unsafe pipeline-contract change

Medium:
- Behavior gap, insufficient validation, blocked required check, or scope risk

Low:
- Maintainability/readability issue or minor inconsistency

For each finding include:

- File path
- Line reference or closest anchor
- Issue summary
- Why it matters
- Recommended minimal fix
- Validation evidence

If no findings:
Write exactly: "No findings."
Then list residual risks and test gaps.
