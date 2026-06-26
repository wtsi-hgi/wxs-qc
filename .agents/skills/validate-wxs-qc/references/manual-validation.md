# Manual Validation Scenarios

## Pipeline Step Changes

1. Identify the numbered stage and its upstream/downstream dependencies.
2. Confirm input paths, output paths, config keys, and Hail object schemas.
3. Run `make integration-test-trios` to test the pipeline end-to-end through the trio integration smoke-test path.
4. To test or debug one stage, run `make test-it-one-step test=test_trios_...` with the concrete trio test name or prefix for that step.
5. If integration smoke tests are blocked, document the missing Hail/Spark/data prerequisites.

## Config Changes

1. Check the public example config and docs for matching keys.
2. Verify defaults and optional fields behave correctly.
3. Confirm path conversion still happens at the Spark/Hail API boundary.

## Shared Helper Changes

1. Search all call sites in numbered scripts, `wes_qc/`, `utils/`, and tests.
2. Confirm no caller relies on the old behavior.
3. Validate with `make integration-test-trios` or the relevant `make test-it-one-step test=test_trios_...` command.
