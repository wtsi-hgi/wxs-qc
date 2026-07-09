---
name: implement-wxs-qc
description: Implement an approved scoped change in WxS-QC, preserving Python/Hail pipeline contracts, numbered stage behavior, config compatibility, and strict scope discipline.
---

# Implement WxS-QC Skill

Use this skill when implementing code after a human has approved a plan.

## Workflow

1. Confirm approval
   - Locate the approved plan in `artifacts/1_plan.md` when one exists.
   - If approval is not explicit, stop and request approval. Do not implement.

2. Re-check repository constraints
   - Read `AGENTS.md`.
   - Confirm the approved scope and files to touch.
   - Re-read the relevant numbered script, helper module, config, and tests before editing.
   - If the approved plan conflicts with the actual codebase, stop and ask for instructions.
   - If required changes fall outside approved scope, stop and ask for confirmation.
   - Assume the environment is already configured and runnable. If tooling, permissions, credentials, data access, cloud access, Spark/Hail setup, or local configuration are missing or broken, stop and report the blocker instead of working around it.

3. Implement minimal changes
   - Keep the diff small, focused, and behavior-preserving unless behavior change is explicitly approved.
   - Preserve numbered script order, public script names, CLI arguments, config keys, and output paths unless explicitly approved.
   - Follow the style of the touched files.
   - Do not apply broad formatting, import cleanup, or refactoring outside the requested code path.

4. Preserve Hail/Spark contracts
   - Keep Hail Table/MatrixTable loading and saving in the same layer unless the approved plan changes that boundary.
   - Prefer passing Hail objects through pipeline functions instead of introducing hidden IO.
   - Treat helper functions as accepting lazy Hail objects by default unless their docstring says otherwise.
   - When a function requires or strongly recommends materialized Hail input, or materializes its output for correctness or performance, document that explicitly in the docstring with `Input contract:` and/or `Output contract:` sections.
     Example: `Input contract: Recommended materialized input MatrixTable. Output contract: Materialized Table.`
   - For normal lazy transforms, either omit special contract text or state the lazy return behavior in `Returns:`, including that the caller owns checkpointing or writing.
   - Convert paths to Spark/Hail format only where the Spark/Hail API requires it.
   - Be careful with keying, row/column annotations, checkpointing, repartitioning, and expensive actions.

5. Update related surfaces only when in scope
   - If config behavior changes, update affected config examples and docs only when approved.
   - If shared helpers in `wes_qc/` or `utils/` change, check all known call sites.
   - Add or update focused tests when the touched area has an existing pattern.

6. Run implementation checks
   - Run checks against modified files:
       - `make check`
       - `make typecheck`
   - Run the individual trio integration test for each changed pipeline step:
     - `make test-it-one-step test=test_trios_<step-name-or-prefix>`
   - Use the concrete `test_trios_...` target that maps to the changed step or the narrowest available step prefix.
   - Do not run the full end-to-end integration suites from the implementer role. The user owns those long-running checks:
     - `make integration-test-trios`
     - `make integration-test-non-trios`
   - If Hail/Spark, cloud data, permissions, missing tools, credentials, or local configuration block checks,
     stop the affected check path and record the blocker clearly.
     Do not modify code, install substitutes, or change configuration to bypass the environment issue.

7. Prepare implementation notes
   - Save notes to `artifacts/2_implement.md` when using the artifacts workflow.
   - Include what changed, why, files touched, checks run, blocked checks, risks, and any out-of-scope discoveries.
   - Read `references/change-checklist.md` before concluding.
