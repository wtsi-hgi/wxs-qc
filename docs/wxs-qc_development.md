# Developer's howto

This howto contains development workflow and best practices for the WxS-QC pipeline.

WxS-QC is an old Python/Hail pipeline with a long history of development and adoption.
It contains numbered pipeline stages, standalone scripts, shared helper modules, Hail/Spark-specific behavior, docs,
and tests that cover only a limited part of the actual behavior.

Treat the existing code as the source of truth,
and perform iterative improvement to align it with coding best practices and style guidance.

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
- `tests/`: unit and integration tests.


## Development environment

### Setting up enviroment for local development

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

Run Ruff and mypy on every modified Python file,
plus the other applicable pre-commit checks on all modified files, through the repository target:

```bash
make check
```

A separate target runs mypy on every modified Python file:

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
```

## Development and code organization best practices

This section contains major suggestions to maintain code style and structure for the WxS-QC pipeline.
These guidelines represent the desired direction for the codebase, but not all pipeline parts follow them yet.
Prefer the local style and contracts of the files you are touching.

### Scripts organization and sequence

- Use “stage” for major pipeline blocks, “step” for individual scripts,
- and “substage” for CLI-selectable blocks within a script.
- Use numbered scripts for pipeline steps.
- Break complex steps into smaller substages using command-line arguments when
  that pattern already fits the step.

### Main function structure

- Keep data loading and saving, especially Hail structures, in the `main()` layer when that matches nearby code.
- Follow the `main()` standard function structure for the step being changed.
- Where possible, avoid moving expensive Hail/Spark IO across function boundaries.

### Pipeline step function design

- Prefer pipeline step functions that accept and return Hail objects.
- Avoid hidden MatrixTable or Table reads and writes inside helper functions.
- If needed checkpoint intermediate results inside functions to a temporary location using `hail.utils.temp_file()`.
- Use dictionary from parsed config for flexible argument passing.
  If needed, unpack the `config` dictionary into individual arguments when it matches the existing call pattern

```python
fstat_hist = plot_f_stat_histogram(sex_ht, **config["stage2"]["f_stat_outliers"])
```

- Convert file paths to Spark format, such as `file://`, only at the point where
  a Hail or Spark API requires that format.

### Toolset organization

- Maintain reusable utility functions in the `wes_qc` package.
- Separate service modules by purpose.

### Hail objects contract documentation

All Hail objects are lazy and do not materialize until they are used.
By default, all helper functions accept lazy Hail objects unless their docstring says otherwise.

However, some functions may need to use materialized Hail objects as input
or materialize outputs to perform their job efficiently.
If a function requires or strongly recomments materialized input for performance or correctness,
it should be annotated with the docstring describing the contract:

"""
  Input contract:
      Recommended materialized input MatrixTable.
      The function scans `mt` multiple times, and assumes the input MatrixTable has already been checkpointed by the caller.
  Output contract:
      Materialized Table.
"""

For functions that are normal lazy transforms, omit it or write only:

  """
  Returns:
      Lazy MatrixTable. The caller owns checkpointing or writing.
  """

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
