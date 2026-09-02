.PHONY: test-unit check typecheck

# Prepend the repository root, so that `wxs_qc`, `utils` and the pipeline stages are importable
export PYTHONPATH:=$(if $(PYTHONPATH),$(PYTHONPATH):)$(shell pwd)
export PYSPARK_PYTHON:=$(shell which python)
export PYSPARK_DRIVER_PYTHON:=$(shell which python)

# Results root written by the integration tests.
# Must stay in sync with general.data_root in tests/integration_tests/config_test_template.yaml,
# i.e. tests/data/<general.dataset_name>_results.
RESULTS_DIR:=tests/data/control_set_small_v2_results

# Runs pre-commit hooks only on modified files
# used to run with agent skills because agents can't stage/commit files
check:
	set -e; \
	tmp_file=$$(mktemp); \
	trap 'rm -f "$$tmp_file"' EXIT; \
	git diff --name-only --diff-filter=ACMR -z HEAD -- > "$$tmp_file"; \
	git ls-files --others --exclude-standard -z >> "$$tmp_file"; \
	if [ -s "$$tmp_file" ]; then \
		xargs -0 pre-commit run --files < "$$tmp_file"; \
	else \
		echo "No modified files to check."; \
	fi

typecheck:
	set -e; \
	tmp_file=$$(mktemp); \
	trap 'rm -f "$$tmp_file"' EXIT; \
	git diff --name-only --diff-filter=ACMR -z HEAD -- '*.py' > "$$tmp_file"; \
	git ls-files --others --exclude-standard -z -- '*.py' >> "$$tmp_file"; \
	if [ -s "$$tmp_file" ]; then \
		xargs -0 mypy --config-file=pyproject.toml < "$$tmp_file"; \
	else \
		echo "No modified Python files to typecheck."; \
	fi

# Runs only the working unit tests; the remaining ones in tests/unit_tests require refactoring
test-unit:
	pytest tests/unit_tests/test_pedigree.py -vv

test-it-one-step:
	cd tests/integration_tests && pytest -vv -s --exitfirst -k $(test)

test-it-one-step-profile: clear-hard-filter-checkpoints
	mkdir -p tests/data
	cd tests/integration_tests && \
	python -m cProfile -o ../data/profile.stats -m pytest -vv -s --exitfirst -k $(test) && \
	snakeviz --hostname 0.0.0.0 $(shell pwd)/tests/data/profile.stats

integration-test-trios: clear-ht clear-logs
	cd tests/integration_tests && pytest -vv -ra --tb=short --exitfirst -k "test_trios_ and not test_trios_0_3_import_data"

integration-test-non-trios: clear-ht clear-logs
	cd tests/integration_tests && pytest -vv -ra --tb=short --exitfirst -k "test_non_trios_ and not test_non_trios_0_3_import_data"


integration-test-coverage: clear-ht clear-logs
	cd tests/integration_tests && pytest -vv -ra --tb=short --exitfirst -k "test_trios_ and not test_trios_0_3_import_data" --cov=../..


clear-hard-filter-checkpoints:
	rm -rf $(RESULTS_DIR)/annotations/testhash/json_dump/* || true

clear-logs:
	rm hail*.log || true
	rm hlrun_* || true
	rm tests/unit_tests/hail*.log || true
	rm tests/integration_tests/hail*.log || true
	rm tests/data/hail*.log || true
	rm $(RESULTS_DIR)/hail*.log || true

clear-ht:
	rm -rf $(RESULTS_DIR)/matrixtables/* || true
