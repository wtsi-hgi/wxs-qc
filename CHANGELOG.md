# Changelog

All notable changes to the WxS-QC pipeline, grouped by version tag.

Versions before `v.0.8` are not covered: the earlier tags in this repository
(`ALSPAC-2023`, `BiB-2023-2024`, `MCS-2023`, `Hardfilter-missingness-ALSPAC-MCS-2024`)
mark dataset-specific analysis snapshots rather than pipeline releases.

---

## [0.9.4] - 2026-09-03

Changes accumulated on `main` after the `v.0.9.3` tag.
The version in `pyproject.toml` is already bumped to `0.9.4`, but the tag is not created yet.

### Pipeline algorithm

* **Variant QC: random test set for the random forest.**
  The RF test set is no longer a fixed genomic interval (`chr20`).
  Instead, a configurable percentage of the labelled TP/FP variants is drawn at random
  (`stage3 -> rf_test_percentage`, default 20%).
  TP and FP variants are sampled separately, so the test set keeps the TP/FP proportions of the full set.
  A new `stage3 -> train_rf -> seed` parameter fixes the random draw for reproducible runs;
  with `null` a new test set is selected on every run.
  The RF run metadata in `rf_runs.json` now records `test_percentage` instead of `test_intervals`.
* **Sample QC: introduced experimental pedigree validation and correction (step 2.2).**
  A new `--validate-trios` substep compares the pedigree against the PC-Relate estimates and
  retains only the trios in which both the maternal and the paternal parent-offspring
  relationships are supported by the kinship and IBD0/IBD1/IBD2 values
  (thresholds defined in the new pedigree threshold constants).
  The corrected pedigree is written to `stage2 -> relatedness_output -> revised_pedigree`
  (set to `null` to disable the correction).
  The revised pedigree can then be used as `stage3 -> pedfile`.
  PC-Relate is now run with `statistics: 'all'` (previously `kin2`), because the correction
  needs IBD0/IBD1/IBD2 in addition to the kinship coefficient.
  The relatedness table is also kept as a Hail Table (`relatedness_output -> relatedness_ht`).

### Bug fixes

* Variant QC: a variant satisfying both the TP and the FP criteria was annotated as
  `tp = true` **and** `fp = true` at once, and so entered RF training in both classes.
  TP/FP labels are now mutually exclusive, and missing truth-set annotations are treated as `false`.
* Genotype QC: fixed the founder-detection logic used by the hard-filter evaluation (step 4.1).
  `select_founders` subtracted *every* proband from the pedigree sample set, so a parent that
  carries its own `.fam` record was wrongly excluded from the founders.
  Founders are now defined as the samples that do not name a parent in the pedigree.
  This corrects the transmitted/untransmitted synonymous singleton ratio (`t_u_ratio`)
  reported by step 4.1, which was computed from a truncated - and sometimes empty - founder set.
  Covered by new unit tests.
* Sample QC: batch annotation no longer depends on sorting samples by name;
  a two-pass annotation is used instead, which produced wrong batch assignments for some datasets.

### Technical improvements

* **New fixed test dataset (`control_set_small_v2`).**
  `scripts/testdata-create/` now contains the full generation chain for a reproducible test dataset:
  CRAM download, Sarek manifest and variant calling, sex and trio metadata, VerifyBamID,
  and VEP annotation. Two GIAB trios were added by slicing exome-target reads out of the
  whole-genome BAMs, and one sample was dropped because its name clashed with a sample in the
  1000 Genomes resource dataset. Integration tests, configs, paths, and the Makefile were
  migrated to the new dataset.
* `select_founders` and `collect_pedigree_samples` were moved out of `utils/utils.py` into a new
  `wxs_qc/pedigree.py` module, with a `make test-unit` target for the working unit tests.
* Obsolete standalone scripts were removed (gtcheck helpers, PLINK/array conversion shell scripts,
  `verifybamid.sh`, `run_tests_remote`), superseded by `scripts/testdata-create/`.
* Documented the `--validate-trios` substep and the `rf_test_percentage` / `seed` parameters,
  described the caveats of relatedness detection in small cohorts,
  and added citation information for the WxS-QC paper.

### Minor

* Added a validation of the `rf_test_percentage` value (must be in `(0, 100]`).
* Added unit tests for the pedigree helper functions.
* Updated the step 2.2 integration test expectations for the pedigree correction.

---

## [v.0.9.3] - 2026-08-04

Genotype QC evaluation improvements and a large technical cleanup:
the package and pipeline stages were renamed, and the test and developer tooling was reworked.

### Pipeline algorithm

* **Genotype QC: validated TP/FP metrics and F1 score in hard-filter evaluation (step 4.1).**
  The evaluation can now take an optional sample-level list of manually validated true- and
  false-positive calls (`general -> metadata -> validated_variants_tsv`, columns
  `type, sample, chr, pos, ref, alt`) and reports, for every filter combination, how many
  validated TP and FP genotypes survive filtering.
  Results are written to `validated_variants_check_tsv` and plotted
  (`snv-TP_validated-FP_validated`, `indel-TP_validated-FP_validated`).
  Precision/recall against the GIAB truth set is now complemented by an F1 score
  (overall and per frameshift/inframe class); `-1` marks combinations where F1 is undefined.
* **Genotype QC: GIAB VCF is now optional (step 4.1).**
  Without a GIAB benchmark VCF the step runs and reports all other metrics,
  skipping only precision/recall.
* **Genotype QC: variant-level filter tags in the output.**
  Step 4.2 annotates variants with a `filters` set
  (`stringent_pass`/`stringent_fail`, `medium_pass`/`medium_fail`, `relaxed_pass`/`relaxed_fail`),
  and the VCF export writes the matching `FILTER` header descriptions,
  so the filtering level of each variant is readable directly from the exported VCFs.
* **Genotype QC: step 4.1 now filters out the variants not used for the metrics calculation**
  before counting, instead of counting on the full matrix table.
* **Sample QC: chrY call rate is computed on male samples only,**
  which removes the systematic call-rate bias caused by including female samples.

### Bug fixes

* Fixed the path formatting for `validated_variants_check_tsv` in the config files.
* Fixed an incorrect variant row count in the filtering log messages.

**Known issue:** the transmitted/untransmitted synonymous singleton ratio (`t_u_ratio`)
reported by step 4.1 is unreliable in this version because of a bug in the pedigree
founder detection. Fixed in 0.9.4.

### Technical improvements

* **Package and stage renaming.**
  The `wes_qc` package is renamed to `wxs_qc`, and the numbered stage directories and scripts
  are renamed to importable module names, e.g.
  `2-sample_qc/2-prune_related_samples.py` → `stage2_sample_qc/sqc2_identify_related_samples.py`,
  `3-variant_qc/3-train_rf.py` → `stage3_variant_qc/vqc3_train_rf_model.py`.
  Top-level config keys are renamed accordingly: `stepN` → `stageN`.
  Dynamic imports in the tests and scripts are replaced with static ones.
* **Partitioning and performance.**
  A global `general -> n_partitions` parameter is added, and the matrix table is repartitioned
  after the heavy filtering steps (variant filtering, sex imputation, step 4.1).
  Repartitioning is done with `shuffle=False` to preserve variant ordering and avoid the
  shuffle overhead; shuffled partitioning was also found to affect Hail LD pruning results.
  `ld_prune` now receives an explicit `memory_per_core` setting.
  Redundant checkpoints and a duplicated autosomal-biallelic-SNV filtering pass were removed,
  and extra checkpoints were added where they speed up repeated filtering.
* **Step 2.2 modularised** into separate substeps (`--filter-mt`, `--pc-relate`, `--plot-pca`, `--all`),
  which makes the step restartable and testable substep by substep.
  `prune_pc_relate` lost its dataset-specific special cases and now enforces schema compatibility.
  Sex imputation was refactored to a cleaner unphased-GT handling with explicit intermediate checkpoints.
* **Test infrastructure.**
  Integration tests migrated from `unittest` to `pytest` with shared fixtures in `conftest.py`;
  the dynamic `__new__`-based test generation was replaced by an explicit ordered step registry.
  Added a `run_integration_remote` script for running the integration suite on remote VMs,
  and reorganised the test data directory layout.
* **Developer tooling.**
  New `make check` (pre-commit on modified files) and `make typecheck` (mypy on modified files)
  targets, expanded pre-commit hooks, ruff and strict mypy configuration moved into `pyproject.toml`.
  Hail updated to 0.2.138 (requirement `>=0.2.137`) to match the cluster setup.
  Added agent role and skill configuration files for AI-assisted development
  (`AGENTS.md`, `.claude/`, `.agents/`, `.codex/`).
* Documented the Hail object input/output contract conventions for pipeline step functions,
  unified the "stage"/"step" terminology across the documentation, added a step 2.2 substage
  description, and explained the content of the exported VCFs in the how-to.

### Minor

* Added output validation for step 2.2 (sample QC) and step 4.1 (genotype QC) in the integration tests.
* Added tests for the alternative sample QC methods (`nn`, `lr`).
* Added a CSV/TSV comparison utility with float tolerance, and an `assert_saved_tables_match` helper;
  the relatedness validation now compares TSV files instead of JSON.
* Added a script to generate the validated TP/FP fixture from an annotated matrix table.
* Added logging of variant counts before and after filtering, of validated TP/FP counts,
  a variant filter statistics printout in the hard-filtering steps, and a warning when
  pedigree founders are unavailable.
* Added a check that the relatedness matrix contains no unnecessary entries before removing samples.
* Standardised the sample field name to `s` across the genotype QC scripts, and added type annotations
  across the utility modules.
* Genotype QC validation is disabled in the step 4.1 integration test because the random forest
  output is non-deterministic.

---

## [v.0.9.2] - 2026-05-26

A configuration and documentation consolidation release on top of `v.0.9`.

### Pipeline algorithm

* Step 2.1: the sex-annotation hard filter threshold `n_alt_alleles_threshold` was raised
  from `0.001` to `0.01`, restricting the variant subset used for sex imputation to more
  common variants.
* Step 4.4 (counts per sample) now returns early with an explicit message when no VEP
  consequences file is provided, instead of silently running the surrounding logic.

### Technical improvements

* The config structure was merged with the sample QC methods introduced in `v.0.9`.
* PCA parameters in the config files were aligned with the `pc_relate_params` constraints:
  `n_pcs` for the nearest-neighbours and linear-regression methods must not exceed
  `pc_components` used for PC-Relate. Documented in the config comments.
* Added the `control_set_small.csq_header.txt` metadata file for the test dataset.
* Documentation: described the new relatedness-check and PCA-based superpopulation prediction
  approaches, the new consequences file structure for VEP annotations and the corresponding
  filtering configuration, the variant filtering parameters used for sample stratification,
  refined the cohort size recommendations in the README, and clarified that the HWE p-value
  is calculated by the Hail `variant_qc` function.

---

## [v.0.9] - 2026-04-28

The main sample QC and genotype QC methods release.

### Pipeline algorithm

* **Three selectable sample QC methods (step 2.4).**
  The outlier detection method is chosen with `general -> sample_qc_method`:
  * `pop` — stratified QC within PCA-assigned superpopulations (the previous behaviour, gnomAD v3);
  * `nn` — nearest neighbours in PC space (gnomAD v4), which handles admixed and
    under-represented ancestries that do not fit the 1000 Genomes superpopulations;
  * `lr` — linear regression residuals of the QC metrics on the principal components (gnomAD v4).
  All methods flag samples deviating from their reference group by more than the configured
  number of median absolute deviations; the MAD table is written out (`mad_file`) and the
  plotting was reworked to plot MAD-based metrics with per-metric thresholds.
  `lr` is the new default in the example config.
* **PC projection for superpopulation prediction (step 2.3).**
  PCA is now computed on the 1000 Genomes reference samples only, and the study samples are
  projected onto those principal components, instead of running PCA on the merged dataset.
  This keeps the PC axes stable and undistorted by related samples or unusual ancestry,
  makes PCs comparable between studies, and removes the limit on the number of related
  individuals in the cohort. Merging and LD pruning were reordered accordingly,
  and the affected output files renamed.
* **KING pre-filtering before PC-Relate (steps 0.1 and 2.2).**
  Related samples are now removed with KING (`kinship_threshold: 0.0442`) before the
  PCA that feeds PC-Relate, both for the 1000 Genomes reference panel and for the study dataset.
  PC-Relate then runs on the merged matrix of related and unrelated samples, so that it uses
  the same variants as the PCA it depends on.
* **Batch-aware sample QC (step 1.2, step 2.4).**
  An optional batch metadata file (`general -> metadata -> batch_file`) can be supplied;
  sample QC statistics are then computed per batch. Samples are marked as `batch1` when
  no batch information is available. Controlled with the `use_batch` parameter.
* **Control samples are excluded from sample QC (steps 2.2-2.5).**
  Samples listed in `general -> metadata -> control_samples` are kept in the dataset
  regardless of the sample QC verdict, so that they remain available for variant QC.
  The final filtering step picks up the matrix table without population annotation.
* **Sex-chromosome-specific genotype hard filters (step 4.2).**
  With `sex_chromosome_specific_filtering: True`, the depth threshold is halved on the
  non-PAR X and Y regions for male samples, using a resolved sex annotation file.
* **VCF export changes (steps 4.3a and 4.3b).**
  Variants with `AC = 0` after filtering are now removed from the exported VCFs;
  step 4.3b can export a VCF for any single specified filter level, not only the stringent one;
  the VEP annotation file and the CSQ annotation of the output VCFs are now optional
  (`cqfile` / `csq_header` may be `null`).

### Bug fixes

* Fixed the crash on step 4.1 when the `assigned_pop` annotation was missing from the matrix table.
* Freemix outlier filtering no longer drops samples with a missing Freemix score.
* Fixed the column index used to read the variant consequence from the VEP consequences file
  (`f5` → `f6`, following the new file structure) and added an explicit filter to
  `synonymous_variant` when selecting synonymous variants.

### Technical improvements

* In the variant QC annotation step, the branching on the presence of the pedigree and of the
  VEP consequences file was moved out of a helper wrapper into `main()`, following the
  convention of keeping IO and workflow branching in the step entry points.
  The trios-related annotation now requires both the pedigree and the consequences file.
* Duplicated relatedness and PCA helper functions were unified into `wes_qc/compute_relatedness.py`
  and `wes_qc/pca_utils.py`; step 2.3 was reworked to use `pca_utils`.
* Relatedness, LD pruning, and PC-Relate parameters were restructured in the config into explicit
  `filter_params`, `king_params`, `prune_params`, `pc_relate_params`, and `relatedness_output` sections,
  with the underlying Hail arguments passed through as keyword dictionaries
  (`ld_prune_args`, `pc_relate_args`).
* Added `utils/mtstat.py`, a small utility to print sample statistics of a matrix table.
* Documentation: substantially extended the how-to with descriptions of the three sample QC
  methods, the relatedness step, the PCA projection approach, and their limitations;
  added the new methods description to the concepts document.

---

## [v.0.8] - 2025-10-14

Baseline release for this changelog. Notable changes in this version:
fixed loading of pre-existing random forest models, and updated Hail and gnomAD versions.
Earlier changes are not tracked here.
