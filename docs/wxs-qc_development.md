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

#### Known gap: relatedness detection is not covered on the small test cohort

Testing relatedness detection was one of the reasons for building a test dataset with trios,
but that coverage is currently disabled and the behaviour is **under validation**.

On `control_set_small_v2` (14 samples, four declared trios) PC-Relate reports kinship around
0.10 for known parent-offspring pairs instead of the expected ~0.25.
Nothing then clears the default `relatedness_threshold` of 0.125, so
`mt_related_samples_to_remove.ht` is empty. Because `run_population_pca` splits the
MatrixTable on that table, every population PCA output covers all 14 samples and the
PC-projection branch for related samples runs on no samples at all.
The approach was validated on large cohorts, so cohort size is the suspected cause,
but this has not been confirmed.

Consequences for validation work:

- In `assert_step_2_2_outputs_match_expected`
  (`tests/integration_tests/test_integration_trios.py`) the assertions downstream of the
  PC-Relate kinship threshold are commented out rather than frozen as correct: the two
  relatedness TSV comparisons, the samples-to-remove count, and all four population PCA
  outputs. The KING-split assertions upstream of the threshold are still active.
- Their expected values are kept in `expected_integration_test_results.json`
  and the reference tables are kept in `tests/integration_tests/validation/`,
  so the checks can be restored unchanged once the behaviour is understood.
- Treat step 2.2 as a smoke test for the related-samples path until then, and do not
  take a passing `test_trios_2_2_sample_qc` as evidence that relatedness detection works.

See the [PC-Relate section of the user howto](wxs-qc_howto.md) for the operational guidance.

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

## Hail Matrixtable I/O operations and optimization

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


## Hail MatrixTable partitioning and performance guide

Hail distributes `Table` and `MatrixTable` data into **partitions**.
Each partition can be processed independently by a Spark task,
allowing Hail to execute computations in parallel across workers and CPU cores.

For large sequencing cohorts, partition count affects
available parallelism, memory required by individual tasks, I/O efficiency, etc.

Hail recommends having approximately **2–4 partitions per available core** as a general rule of thumb.

### How a MatrixTable is partitioned

Consider the usual genetics representation:

```text
                    samples
              ──────────────────>

variant 1     GT  GT  GT  GT  ...
variant 2     GT  GT  GT  GT  ...
variant 3     GT  GT  GT  GT  ...
...
```

A Hail `MatrixTable` is partitioned **only along the row (variant) axis**,
where partitions contain all samples for a contiguous, non-overlapping range of rows.
Hail does **not** normally partition the matrix along sample axis.

Therefore:
* The number of variants controls the partitioning axis.
* The number of samples controls how wide and expensive each row is.
  Samples affect partition size, but not partition parallelism.


### Row keys determine partition ranges

Genomic MatrixTables are normally keyed by something like:

```python
(locus, alleles)
```

Partitions therefore correspond to ranges of genome sorted by varuant locus.

The partitions do not necessarily contain exactly the same number of variants,
because variant density differs across genomic regions.

### Entry schema also affects partition cost

Two MatrixTables with identical dimensions can have very different computational costs.

For example, entry containing only `GT` is much smaller than entry with
`GT,GQ,DP,AD,PL`

Fixed-width fields such as `GT` and `GQ` have relatively predictable costs.

Variable-length fields such as `AD`
can make individual rows significantly wider, especially for multiallelic variants.

A useful conceptual approximation is therefore:

$$
\text{partition workload}
\approx
N_{\text{rows per partition}}
\times
N_{\text{samples}}
\times
B_{\text{average entry}}
$$

where \(B_{\text{average entry}}\) is the approximate amount of entry data per sample and variant.

### Partitions and Spark parallelism

Normally, one Spark task processes one partition at a time.
Suppose a cluster has 50 workers and 8 cores / worker,
giving approximately 400 available cores,
and the MatrixTable has 12,000 partitions.
Then approximately 400 partitions can be processed simultaneously.
Having more partitions than cores is generally desirable because it:
keeps cores occupied, compensates for partitions that take longer than others,
and reduces the impact of stragglers.
Hail currently recommends approximately **2–4 partitions per core** as a basic parallelism guideline.

If a Matrixtable has too few partitions,
a part of the cluster cannot contribute to the processing of this matrixtable.
Increasing partition count can help in this situation.

The opposite extreme is also inefficient.
Small partitions may contain very little useful work, whereas
Spark still has to pay per-task costs for scheduling, task startup, serialization, etc.

There is no universally optimal partition count
because different operations have different memory and compute requirements.
A useful starting heuristic is:

$$
P_{\text{target}}
\approx
\max
\left(
2\text{–}4\times N_{\text{cores}},
P_{\text{size}}
\right)
$$

where \(P_{\text{size}}\) is the number required to keep individual partitions at a reasonable data size.

For pipeline heuristics, a practical starting target is often around ~128 MiB logical data per partition
This size target is a **practical heuristic**, not a strict Hail requirement. It should be adjusted according to
entry width, operation type, memory per core, etc.

### Filtering changes what the existing partitions mean

If you perform a heavy filtering and remove about significant part of variants,
the partition number of MatrixTable remains the same,
although each partition now contains much less data (or even be empty).
Hail explicitly notes that partition overhead means reducing partition count can improve performance
after significant filtering.

Therefore, you usually should reduce the number of partitions,
depending on the number of filtered variants.

Keep in mind that row filtering and column filtering behave differently.
Row filtering removes variants.
Because rows define partitions, heavy row filtering can leave many undersized or empty partitions.
Column filtering removes samples.
This makes every MatrixTable row narrower but does not divide or remove row partitions in the same way.
Therefore major sample filtering can also justify reducing partition count,
even though the number of rows has not changed.
The same logic applies to dropping entries, because it also makes partitions narrower.


### Repartitioning

Repartitioning with `shuffle=False`
combines existing partitions and avoids a full data shuffle.
Hail documents this as analogous to Spark `coalesce`.
Consequently, it can be applied only if
`target_partitions <= mt.n_partitions()`
This is usually appropriate when existing partitions are reasonably balanced,
and filtering is distributed throughout the genome. For example filtering for
biallelic variants, SNVs, PASS variants, and common genome-wide QC filters.

Calling `repartition(..., shuffle=True)`
performs a full shuffle and creates more evenly sized partitions.
This is considerably more expensive, so use it when out need to
increase partition count, fix with badly imbalanced partitions or
highly irregular genomic regions.
For example, filtering the whole genome matrixtable to only one chromosome,
or small collection of exome intervals.

### VDS and MatrixTable partitioning

A Hail VDS contains two MatrixTables
(reference_data and variant_data), that use the same row-oriented partitioning model.
After `mt = hl.vds.to_dense_mt(vds)`
the resulting dense MatrixTable is also partitioned along variants.

This is important for large WGS cohorts because converting a sparse VDS into a dense MT
can greatly increase the amount of data associated with each row while preserving row-oriented parallelism.
Whenever possible, large variant reductions should therefore be performed on the VDS before densification.

### Different operations can prefer different partition sizes

There is no single partition count that is optimal for every Hail operation.
For example, Row-only operations, like `mt.count_rows()`
do relatively little work per row.

Entry-heavy aggregations, like `hl.agg.mean(mt.DP)`
may inspect millions of entries per partition.
Smaller partitions can help in this case by reducing memory pressure,
improving load balancing, and creating more parallel tasks.

Shuffle-heavy operations, like joins, grouping, repartitioning, and some aggregations
may generate substantial intermediate data.
Their optimal partition count may differ from that of straightforward scans.

### Inspect Spark UI when tuning

Partition heuristics are useful, but Spark UI provides the best evidence about an actual stage.
Having ling tasks (many minutes) and a smaller number of tasks than available cores,
likely indicated problems with insufficient parallelism and oversized partitions.

Opposite, if you see tasks that are complete in a few seconds,
you likely have too small partitions and significant scheduling overhead.

If most tasks finish quickly, but a few tasks take dramatically longer
you probably have uneven partition sizes or unusually expensive genomic regions.


## Recommended pipeline practices

### Record partition metadata

At important pipeline boundaries, record:

```python
n_rows = mt.count_rows()
n_cols = mt.count_cols()
n_partitions = mt.n_partitions()
```

Prefer to calculate counts once and store them in pipeline metadata
rather than repeatedly triggering Hail actions.


### Reconsider partitions after major transformations

Check and change the number of partitions after:

* large variant filtering;
* large sample filtering;
* selecting/dropping large entry fields;
* restricting to chromosomes or intervals;
* joining or grouping datasets;
* producing a substantially smaller QC dataset.

Small transformations normally do not justify repartitioning.


| Situation                                 | Recommended action                       |
| ----------------------------------------- | ---------------------------------------- |
| `<2 partitions/core`                      | Consider increasing partitions           |
| `2–4 partitions/core`                     | Good minimum parallelism                 |
| Large partitions / OOM                    | Increase partitions                      |
| Thousands of tiny fast tasks              | Reduce partitions                        |
| Heavy row filtering                       | Recalculate/reduce partitions            |
| Heavy column filtering                    | Re-estimate because rows became narrower |
| Uniform genome-wide filtering             | Prefer `shuffle=False` when reducing     |
| Highly uneven genomic filtering           | Consider `shuffle=True`                  |
| Need more partitions than currently exist | Requires `shuffle=True`                  |
| Current count is already near estimate    | Do nothing                               |
| Before expensive repartition              | Benchmark/inspect Spark UI               |


### Avoid chasing an exact optimum

Subtle changing of partition number, like 6000 → 5000
is unlikely to pay for itself.
Repartition when the difference is substantial.

## Efficient Variant Filtering and Repartitioning in Hail

Excessive partitions increase task scheduling and I/O overhead,
so it is often useful to reduce the partition count after heavy filtering.
However, a naïve workflow such as:

```python
filtered = mt.filter_rows(...)
n_rows = filtered.count_rows()

filtered = filtered.repartition(...)
filtered.write(...)
```

can evaluate the filtering pipeline twice:
once for `count_rows()` and then again for `write()`.

The preferred strategy is, therefore:

1. Determine which **row keys** should survive.
2. Store those keys in a small `Table` called, for example, `keep_rows`.
3. Count `keep_rows` to estimate the appropriate final partition count.
4. Apply the same `keep_rows` table to the full `MatrixTable`.
   Hail provides `MatrixTable.semi_join_rows()` specifically for keeping MatrixTable rows
   whose keys occur in another keyed `Table`.
5. Repartition/coalesce and write the result.

The optimal implementation differs
depending on whether filtering uses only row information or requires aggregation over entries.

### Filtering based only on row data

Examples include filtering by `locus`, alleles, chromosome,
SNP/indel type, number of alleles, VCF `FILTER`, etc.
This condition does not require examining `GT`, `GQ`, `DP`, or any other entry field.

In this case, instead of filtering the complete MatrixTable immediately, operate on its row table:

```python
rows = mt.rows()

keep_rows = (
    rows
    .filter(
        (hl.len(rows.alleles) == 2)
        & hl.is_snp(rows.alleles[0], rows.alleles[1])
    )
    .select()
)
```

`select()` with no arguments drops non-key fields while retaining the table key.
For a conventional variant MatrixTable this produces a small table conceptually containing only locus and alleles.
rather than genotype data for thousands of samples.

Because `keep_rows` contains only row information,
evaluating it is much cheaper than processing all MatrixTable entries.
If `keep_rows` will be counted and then reused for filtering, you can materialize it.

Now the exact number of retained MatrixTable rows is simply:

```python
n_after = keep_rows.count()
```

If the original number of rows is not known, obtain it:

```python
n_before = mt.count_rows()
```

This is still a row operation and does not require an aggregation over all genotype entries.

Next, estimate the new partition count.
Since filtering changes only the number of variants
while leaving the entry schema approximately unchanged, a useful first approximation is:

\[
P_{\mathrm{new}}
\approx
P_{\mathrm{old}}
\times
\frac{N_{\mathrm{after}}}{N_{\mathrm{before}}}
\]


The simple scaling assumption works best when the average entry size stays approximately constant.
It can be less accurate for fields such as AD or LGT,
because filtering multiallelic variants to biallelic variants may also substantially reduce average array lengths.
In that situation, use row-count scaling as a first estimate or sample retained rows when estimating entry width.

The result should also respect cluster parallelism. Hail recommends approximately 2–4 partitions per available core.

Partition selection is inherently approximate:
the difference between 4,500 and 5,000 partitions is generally not worth an additional complete MatrixTable scan.


Next, apply the precomputed `keep_rows` to the MatrixTable:

```python
filtered_mt = mt.semi_join_rows(keep_rows)
```

`semi_join_rows` keeps MatrixTable rows whose row key appears in the supplied table.
If `keep_rows` was derived from `mt.rows()`, the keys are already compatible.

An equivalent explicit lookup is possible:

```python
filtered_mt = mt.filter_rows(
    hl.is_defined(
        keep_rows.index(mt.locus, mt.alleles)
    )
)
```

but `semi_join_rows()` expresses the intention more clearly.

Finally, repartition and write

After a substantial row reduction, the target partition count will normally be **lower** than the existing partition count.
Therefore, you can repartition without shuffle.

```python
filtered_mt = filtered_mt.repartition(
    target_partitions,
    shuffle=False,
)

filtered_mt.write(
    output_path,
    overwrite=True,
)
```

With `shuffle=False`, Hail combines existing partitions.
This is particularly suitable for filters distributed relatively uniformly across the genome, such as:

```text
biallelic variants
SNVs
PASS variants
common QC categories
```

### Filtering that requires entry aggregation

Some filters cannot be determined from row annotations alone.
Examples include call rate, allele frequency calculated from the cohort, mean or median DP, etc.

In this case, some expensive computation is unavoidable.

#### Aggregation strategy A: materialize only `keep_rows`

This is often the best approach when the aggregated statistics are needed only to decide whether a variant survives.
Use this strategy when:

- QC statistics are needed only for filtering;
- the original native MatrixTable is relatively cheap to reread;
- the filtered MatrixTable is still large;
- minimizing temporary disk usage is important.

Compute the required row metrics and immediately take rows:

```python
row_metrics = (
    mt
    .select_rows(
        call_rate=hl.agg.fraction(
            hl.is_defined(mt.GT)
        ),
        mean_dp=hl.agg.mean(mt.DP),
        mean_gq=hl.agg.mean(mt.GQ),
    )
    .rows()
)
```

Hail permits column/entry aggregation when constructing MatrixTable row fields.

Filter those metrics:

```python
keep_rows = (
    row_metrics
    .filter(
        (row_metrics.call_rate >= 0.99)
        & (row_metrics.mean_dp >= 20)
        & (row_metrics.mean_gq >= 30)
    )
    .select()
)
```

Then checkpoint the **small key table**:

```python
keep_rows = keep_rows.checkpoint(
    "tmp/qc_keep_rows.ht",
    overwrite=True,
)
```

This action performs the expensive aggregation.
Subsequent operations use the stored result rather than recomputing the aggregation.

Now count rows in the `keep_rows` table and calculate the desired partition count.

Finally, run `semi_join`, repartition, and save the result:

```python
filtered_mt = mt.semi_join_rows(keep_rows)

filtered_mt = filtered_mt.repartition(
    target_partitions,
    shuffle=False,
)

filtered_mt.write(
    output_path,
    overwrite=True,
)
```

#### Aggregation strategy B: checkpoint the filtered MatrixTable

Sometimes the computation upstream of filtering is sufficiently expensive that rereading/reconstructing the MatrixTable is undesirable.

In that case, materialize the filtered MT itself:

```python
filtered_mt = (
    mt
    .annotate_rows(
        call_rate=hl.agg.fraction(
            hl.is_defined(mt.GT)
        ),
        mean_dp=hl.agg.mean(mt.DP),
    )
    .filter_rows(
        lambda row: ...
    )
)

filtered_mt = filtered_mt.checkpoint(
    "tmp/filtered.mt",
    overwrite=True,
)
```

`MatrixTable.checkpoint()` materializes the MatrixTable using a fast, less space-efficient codec and rereads it.
It is intended to break long computation graphs and avoid repeated upstream computation.

Now obtain the exact row count from the checkpoint,
choose the desired partition count, and write the optimized final result.

Use this strategy when:

- the filtering aggregation is very expensive;
- there is a long or expensive upstream computation graph;
- the aggregate row fields should remain in the final MatrixTable;
- rereading/recomputing the original source would be expensive;
- temporary disk space is available.


#### Keeping QC metrics without checkpointing the whole MT

There is also a useful hybrid approach. Instead of retaining only variant keys,
checkpoint the complete `row_metrics` table.
Then construct `keep_rows` cheaply from that checkpoint:

```python
keep_rows = (
    row_metrics
    .filter(
        (row_metrics.call_rate >= 0.99)
        & (row_metrics.mean_dp >= 20)
    )
    .select()
)
```

The expensive entry aggregation has already been materialized.
You can then count rows and run `semi_join_rows`

This pattern is especially useful when several filtering thresholds will be tested against the same expensive QC metrics.
