.PHONY: test results all

test:
	python -m unittest discover -s tests -v

results:
	python -m numerics.run_compare --points 120

all: results test
