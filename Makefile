.PHONY: test check typecheck

export PYTHONPATH:=$PYTHONPATH:$(shell pwd)
export PYSPARK_PYTHON:=$(shell which python)
export PYSPARK_DRIVER_PYTHON:=$(shell which python)

# Runs pre-commit hooks only on modified files
# used to run with agent skills because agents can't stage/commit files
check:
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
	tmp_file=$$(mktemp); \
	trap 'rm -f "$$tmp_file"' EXIT; \
	git diff --name-only --diff-filter=ACMR -z HEAD -- '*.py' > "$$tmp_file"; \
	git ls-files --others --exclude-standard -z -- '*.py' >> "$$tmp_file"; \
	if [ -s "$$tmp_file" ]; then \
		scripts/stage_mypy_numbered_scripts.sh; \
		xargs -0 mypy --config-file=pyproject.toml < "$$tmp_file"; \
	else \
		echo "No modified Python files to typecheck."; \
	fi

test-it-one-step:
	cd tests/integration_tests && pytest -vv -s --exitfirst -k $(test)

test-it-one-step-profile: clear-hard-filter-checkpoints
	cd tests/integration_tests && \
	python -m cProfile -o profile.stats -m pytest -vv -s --exitfirst -k $(test) && \
	snakeviz --hostname 0.0.0.0 $(shell pwd)/tests/integration_tests/profile.stats

integration-test-trios: clear-ht clear-logs
	cd tests/integration_tests && pytest -vv -ra --tb=short --exitfirst -k "test_trios_ and not test_trios_0_3_import_data"

integration-test-non-trios: clear-ht clear-logs
	cd tests/integration_tests && pytest -vv -ra --tb=short --exitfirst -k "test_non_trios_ and not test_non_trios_0_3_import_data"


integration-test-coverage: clear-ht clear-logs
	cd tests/integration_tests && pytest -vv -ra --tb=short --exitfirst -k "test_trios_ and not test_trios_0_3_import_data" --cov=../..


clear-hard-filter-checkpoints:
	rm -rf tests/integration_tests/integration-data/annotations/testhash/json_dump/* || true

clear-logs:
	rm hail*.log || true
	rm hlrun_* || true
	rm tests/unit_tests/hail*.log || true
	rm tests/integration_tests/hail*.log || true

clear-ht:
	rm -rf tests/integration_tests/matrixtables/* || true
