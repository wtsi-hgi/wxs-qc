# Developer's howto

This howto contains development workflow and best practices for the WxS-QC pipeline.

WxS-QC is an old Python/Hail pipeline with a long history of development and adoption
For now, the pipeline code has limited automated coverage, and workflow contracts that are not always fully documented.
Treat the executable code as the source of truth
and perform iterative improvement to align it with coding best practices and style guidance.

## Development environment

Update your environment with all dependencies, including dev and test packages:

```bash
uv sync
```

Set up `pre-commit` for local commits:

```bash
pre-commit install
```

This installs hooks for checks such as trailing whitespace removal, end-of-file fixes,
YAML syntax checking, large file checking, Ruff linting and formatting, and MyPy type checking.

### Agent-assisted development

The repository contains local instructions for AI-assisted work:

- `AGENTS.md`: canonical repository rules and scope constraints.
- `.agents/roles`: planning implementation and validation roles.
- `.agents/skills`: Corresponding skills for each role.

For agent-assisted changes use the plan-implement-validate approach with a himan in the loop,
which requires explicit approval after the planning step.
Agent implements one change at a time.
Then a developer should manually review the change, run tests is necessary,
and manually commit changes.

### Checks and smoke tests

All available checks are implemented as Makefile targets.

Pre-commit checks through the repository target:

```bash
make check
```

A separate target for type checking:

```bash
make typecheck
```

Run a targeted trio integration smoke test for each changed pipeline step:

```bash
make test-it-one-step test=test_trios_<step-name-or-prefix>
```

Examples:

```bash
make test-it-one-step test=test_trios_1
make test-it-one-step test=test_trios_2_2
```

Unit tests are currently broken and are not the active validation path.

Full end-to-end integration suites are long-running checks:

```bash
make integration-test-trios
make integration-test-non-trios
```

### Pre-commit hooks

Once `pre-commit` is installed, hooks run automatically on commit.

For repository validation, prefer `make check`;
it is the supported entry point for pre-commit checks in the agent workflow and avoids sandbox-specific command issues.

For manual local development outside the agent workflow, direct pre-commit
commands can still be useful:

```bash
pre-commit run --all-files
pre-commit run --files <file1> <file2>
pre-commit run mypy --hook-stage manual
```

## Development and code organization best practices

This section contains major suggestions to maintain code style and structure for the WxS-QC pipeline.
These guidelines represent the desired direction for the codebase, but not all pipeline parts follow them yet.
Prefer the local style and contracts of the files you are touching.

### Scripts organization and sequence

- Use numbered scripts for pipeline steps.
- Break complex steps into smaller substeps using command-line arguments when
  that pattern already fits the stage.

### Main function structure

- Keep data loading and saving, especially Hail structures, in the `main()` layer when that matches nearby code.
- Follow the `main()` standard function structure for the stage being changed.
- Where possible, avoid moving expensive Hail/Spark IO across function boundaries.

### Pipeline step function design

- Prefer pipeline step functions that accept and return Hail objects.
- Avoid hidden MatrixTable or Table reads and writes inside helper functions.
- If needed checkpoint intermediate results inside functions to a temporary location using `hail.utils.temp_file()`.
- Use dictionary from parsed config for flexible argument passing.
  If needed, unpack the `config` dictionary into individual arguments when it matches the existing call pattern

```python
fstat_hist = plot_f_stat_histogram(sex_ht, **config["step2"]["f_stat_outliers"])
```

- Convert file paths to Spark format, such as `file://`, only at the point where
  a Hail or Spark API requires that format.

### Toolset organization

- Maintain reusable utility functions in the `wes_qc` package.
- Separate service modules by purpose.

## Performance optimization tips

Hail and Spark utilize lazy computational approach,
and most transformations are just a recipe until an action forces execution.

Checkpoints are a mechanism to save and reuse expensive intermediate results.
However, writing checkpoints is an expensive action and can be a bottleneck if used without a reason.

**Avoid checkpoints when**:

- You perform a linear set of filtering, aggregation, or join operations.
- You need to calcuate several aggregations over the same data.
  Instead, put all aggregations in a single aggregation call call.

**Checkpoint when** the recipe becomes expensive, reused, unstable, or too large.
Common checkpoint points include:

- Before branching the pipeline into several downstream outputs.
  Example: a filtered MatrixTable used for KING, PC-Relate, PCA, plotting, and exports.
  Without checkpoint/readback, each action may recompute the filter chain.
- After a large shuffle or expensive aggregation (Which may hide inside Hail functions):
  such as `variant_qc`, `sample_qc`, `ld_prune`, `pc_relate`,
  group-by, interval filters over large tables, or repartitions.
- When lineage is getting long. Many chained filters, annotations, joins, semi-joins, unions,
  and aggregations can make the Spark plan huge and fragile.
