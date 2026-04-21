# Makefile for the WSU data-properties pipeline.
#
# Tests and pipeline runs produce a lot of intermediate artifacts
# (realization ecsv files, generated LaTeX tables, PNG plots, JSON
# stats).  This Makefile separates throwaway test runs from durable
# "validation" runs and makes cleanup trivial.
#
# Common flows:
#   make validation         -- reproduce canonical tables + plots in build/validation/
#   make test-keep          -- run pytest; leave artifacts in build/test_artifacts/
#   make test               -- run pytest (auto-cleaned)
#   make clean              -- remove everything in build/

VENV     := /Users/briantest/Code/alma/venv
PY       := $(VENV)/bin/python
PYTEST   := $(VENV)/bin/pytest
DB       := data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv

BUILD       := build
ARTIFACTS   := $(BUILD)/test_artifacts
VALIDATION  := $(BUILD)/validation

# Realization count for validation runs (memo tables used 10).
N_REAL ?= 10

.PHONY: help test test-keep test-fast validation \
        clean clean-test clean-validation clean-realizations list-artifacts

help:
	@echo "Targets:"
	@echo "  test              -- run full pytest (~6 min; artifacts auto-cleaned)"
	@echo "  test-fast         -- run only fast tests (formulas + sample DB + parser)"
	@echo "  test-keep         -- run pytest and preserve artifacts in $(ARTIFACTS)/"
	@echo "  validation        -- run pipeline CLI to produce canonical outputs in $(VALIDATION)/"
	@echo "                       (tables, plots, stats.json, realizations)"
	@echo "  list-artifacts    -- show what is currently in $(BUILD)/"
	@echo "  clean-test        -- remove $(ARTIFACTS)/"
	@echo "  clean-validation  -- remove $(VALIDATION)/"
	@echo "  clean-realizations -- remove data/sample_band1_band2/ (regenerate takes minutes)"
	@echo "  clean             -- remove $(BUILD)/ entirely"
	@echo ""
	@echo "Variables:"
	@echo "  N_REAL=<n>        -- number of realizations for validation (default 10)"

test:
	$(PYTEST) test_pipeline.py -v

test-fast:
	$(PYTEST) test_pipeline.py -v -k "Formula or SampleDB or TableParser"

test-keep:
	mkdir -p $(ARTIFACTS)
	$(PYTEST) test_pipeline.py -v --basetemp=$(ARTIFACTS)
	@echo ""
	@echo "Test artifacts preserved under: $(ARTIFACTS)/"
	@echo "Tip: 'make list-artifacts' to see the layout."

validation: $(VALIDATION)/stats.json
	@echo ""
	@echo "Validation products written to $(VALIDATION)/"
	@echo "  tables/       -- memo-style + SDD-style LaTeX tables"
	@echo "  plots/        -- SDD-style CCDFs (productsize, datavol, datarate, sysperf)"
	@echo "  realizations/ -- per-MOUS ecsvs with Band 1/2 synthesised"
	@echo "  stats.json    -- computed medians / TWA / max / totals per milestone+array"

$(VALIDATION)/stats.json: $(DB) pipeline.py tables.py plots.py config.py
	mkdir -p $(VALIDATION)
	$(PY) pipeline.py \
	    --db-path $(DB) \
	    --outdir $(VALIDATION)/realizations \
	    --n-realizations $(N_REAL) \
	    --stats-out $(VALIDATION)/stats.json \
	    --generate-tables $(VALIDATION)/tables \
	    --generate-plots $(VALIDATION)/plots \
	    -v

list-artifacts:
	@if [ -d $(BUILD) ]; then \
	    find $(BUILD) -maxdepth 4 -type f \( -name '*.tex' -o -name '*.png' -o -name '*.json' -o -name '*.ecsv' \) \
	        | sort; \
	else \
	    echo "(no $(BUILD)/ directory; nothing to show)"; \
	fi

clean-test:
	rm -rf $(ARTIFACTS)

clean-validation:
	rm -rf $(VALIDATION)

clean-realizations:
	rm -rf data/sample_band1_band2

clean:
	rm -rf $(BUILD)
