---
name: plan-change
description: Create a small, reviewable implementation plan for WxS-QC, an old Python/Hail QC pipeline with limited tests, using references/plan-template.md, scoped files, acceptance criteria, and validation steps. Ends by requesting human approval before implementation.
---

# Plan Change Skill

Use this skill when the user asks to plan a feature, refactor, fix, or enhancement.

## Workflow

1. Read the user request and restate the target behavior.
2. Read `AGENTS.md` and enforce all mandatory constraints.
3. Produce a minimal one-PR plan using the section structure from `references/plan-template.md`.
4. Save the complete proposed plan to `artifacts/1_plan.md`.
5. Keep scope tight:
   - Explicitly list out-of-scope items.
   - Avoid broad refactors if not requested explicitly.
6. Mention pipeline-contract impact explicitly:
   - Either "No pipeline contract changes".
   - Or list required coordinated updates to scripts, config keys, docs, tests, and output paths.
7. Include validation commands and manual checks with clear pass criteria.
8. If requested behavior conflicts with the actual codebase, stop and ask for clarification instead of planning around the conflict.
9. Assume required tools, permissions, credentials, data access, cloud access, Spark/Hail setup, and local configuration are already available.
   If inspection reveals an environment blocker, stop and report it instead of planning a workaround.
10. Stop and request human approval. Do not implement.

## Guardrails

- Prefer incremental changes to broad rewrites.
- Prioritize the lowest regression-risk path when uncertain.
- Keep plans focused enough for a practical review.
- Treat missing/blocked validation steps as risks to call out in the plan.
- Do not plan edits outside the approved scope without marking them as requiring user confirmation.
- Do not propose project-code or configuration changes to compensate for missing tools, permissions, credentials, data, cloud access, Spark/Hail setup, or broken local configuration unless explicitly requested.

## Output

Return ONLY the completed plan sections in the exact order and headings from `references/plan-template.md`.
At the very end, add a single line:
Approval needed: Please confirm this plan is approved. I will not implement until you approve.

Also persist the same content to `artifacts/1_plan.md` before finishing.
