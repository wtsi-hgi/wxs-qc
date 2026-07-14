---
name: validate-wxs-qc
description: Validate WxS-QC changes by running integration smoke tests, reviewing pipeline contracts, checking scope discipline, and documenting blocked integration/Hail/Spark validation.
---

# Validate WxS-QC Skill

Use this skill when asked to review or validate implementation quality.

## Workflow

1. Read `AGENTS.md`.
2. Read `artifacts/1_plan.md` and `artifacts/2_implement.md` when present.
3. Review changed files against the approved scope and actual codebase behavior.
4. Validate whether the implemented changes achieve the goal listed in the `artifacts/1_plan.md`. Especially check that:
   - The code logic is consistent and achieved the required goal
   - The touched tests actually test the required behavior
5. If there are any logical issues, report them and stop the validation.
6. Search all call sites in numbered scripts, `wes_qc/`, `utils/`, and tests.
   If there are any parts that rely on the old behavior, report them and stop the validation.
7. Run `.agents/skills/validate-wxs-qc/scripts/validate.sh` and ensure that there ar no validation errors.
8. Run `make test-it-one-step test=test_trios_...` with the concrete trio test name or prefix for each changed step.
9. Save the validation report to `artifacts/3_validate.md`.
10. Report findings by severity with file/line references.
11. If no findings, state that explicitly and note residual risks.
12. Assume the environment is pre-configured and runnable. If permissions, tools, credentials, data access, cloud access, Spark/Hail setup, or local configuration block validation, stop the affected validation path and report the blocker. Do not install substitutes, alter project code or configuration, relax validation, or widen scope to work around the environment issue.

## Required Checks

Use linters and per-step integration smoke tests as the validation baseline:

- `make check`
- `make typecheck`
- `make test-it-one-step test=test_trios_...` for each changed pipeline step, using the same concrete per-step trio tests run by the implementer

**Notes:**

- Unit tests are currently broken and not used for validation.
- Integration tests may require Hail/Spark setup and generated intermediate files.
- Do not run full end-to-end integration suites from the validator role. The user owns these long-running checks and may run them on a separate machine:
  - `make integration-test-trios`
  - `make integration-test-non-trios`
- Treat blocked checks as residual risk unless equivalent evidence exists.
- Passing integration smoke tests do not prove broad pipeline safety because coverage is limited.
- Missing tools, permissions, credentials, data, cloud access, Spark/Hail setup, or broken local configuration are environment blockers to report, not reasons to change repository code or validation scope.

## Review Focus

- Scope matches the approved plan and `AGENTS.md`.
- the implemented changes achieve the goal listed in the `artifacts/1_plan.md`.
   - The code logic is consistent and achieved the required goal
   - The touched tests actually test the required behavior
- No unrelated files, stages, docs, config, or formatting were changed.
- Public scripts, CLI flags, config keys, and output paths remain compatible unless explicitly approved.
- Hail Table/MatrixTable schemas, keys, annotations, and IO behavior are preserved or intentionally changed.
- Shared helper changes have known call sites reviewed.

## Reporting Format

Start with:

Status: PASS / FAIL / BLOCKED

Then list findings ordered by severity:

High:
- Logical issues, regression, data loss, invalid QC result, or unsafe pipeline-contract change

Medium:
- Insufficient validation, blocked required check, mypy type inconsistency, or scope risk

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
